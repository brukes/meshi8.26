/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.optimizers;

import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.util.Terminator;
import meshi.util.UpdateableException;
import meshi.util.Utils;

/**
 * Minimize energy according to a given setResidue of coordinates and an energy function
 */

public abstract class Minimizer extends Optimizer {
    public final int MAX_KICKSTARTS = 10;
    public final double tolerance;
    private double forceMagnitude;
    private int numberOfKickStrarts;
    public static final Terminator terminator = new Terminator();

    public Minimizer(TotalEnergy energy, int maxSteps, int reportEvery, double tolerance) throws UpdateableException {
        super(energy, maxSteps, reportEvery);
        this.tolerance = tolerance;
    }

    public OptimizerStatus run() throws OptimizerException, UpdateableException, EvaluationException {
        return run(true);
    }

    public OptimizerStatus run(boolean testFlag) throws OptimizerException, UpdateableException, EvaluationException {
        OptimizerStatus os = init();
        if (os == OptimizerStatus.CONVERGED) {
            System.out.println("*************** Gradient is to low to minimize *******");
            double e = energy.evaluate();
            double grad = energy().getGradMagnitude();
            System.out.println("Energy is "+e);
            System.out.println("Gradient  is "+grad);
            return OptimizerStatus.CONVERGED;
        }
        numberOfKickStrarts = 0;
        int step;
        for (step = 1; status(step) == OptimizerStatus.RUNNING; step++) {
            boolean minimizationStepOK = minimizationStep();

            if (!minimizationStepOK) {
                if (numberOfKickStrarts >= MAX_KICKSTARTS)
                    throw new OptimizerException("\n\nThe simulation was restarted for " + MAX_KICKSTARTS + " times " +
                            "which is more than allowed.\n" +
                            "So many restarts are indicative of an ill-shaped energy function or " +
                            "an energy differentiation\n");
                try {
                    kickStart();
                    System.out.println("kickstart # " + numberOfKickStrarts + " done");
                }
                catch (OptimizerException oe) {
                    if ((testFlag) && Utils.verbose()) energy.test();
                    throw oe;
                }
                numberOfKickStrarts++;
            }

            if (step % reportEvery == 0)
                Utils.println(energy().report(step));
        }
        return status(step);
    }

    public OptimizerStatus status(int step) {
        if (terminator.dead()) {
            return OptimizerStatus.KILLED;
        }
        forceMagnitude = energy.getGradMagnitude();

        if (forceMagnitude < tolerance) return OptimizerStatus.CONVERGED;
        if (step <= maxSteps) return OptimizerStatus.RUNNING;
        return OptimizerStatus.UNCONVERGED;
    }

    protected abstract OptimizerStatus init() throws OptimizerException, UpdateableException, EvaluationException;

    protected abstract boolean minimizationStep() throws OptimizerException, UpdateableException, EvaluationException;

    protected abstract OptimizerStatus kickStart() throws OptimizerException, UpdateableException, EvaluationException;
}
