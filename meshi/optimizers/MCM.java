/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.optimizers;

import meshi.energy.AbstractEnergy;
import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.energy.pairwiseNonBondedTerms.ExcludedVol.ExcludedVol;
import meshi.energy.simpleEnergyTerms.tether.TetherEnergy;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.scoringFunctions.CombinedEnergyScore;
import meshi.scoringFunctions.MissingFieldException;
import meshi.scoringFunctions.Score;
import meshi.sequences.AlignmentException;
import meshi.util.Logger;
import meshi.util.MeshiProgram;
import meshi.util.UpdateableException;
import meshi.util.Utils;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

public class MCM extends Optimizer {
    protected TotalEnergy energy;
     Minimizer minimizer;
    private TemperatureGenerator temperatureGenerator;
    protected Perturbation perturbation;
    protected Logger log = null;
    protected TetherEnergy allAtomTether = null;
    protected ArrayList<Score> scoreFunctions;
    protected Score optimizationScore;

    public  enum mcmStepResult {
        FAILED, SUCCEEDED, AFTER_PERTURBATION;
        private  int lastSuccess;
        protected static void setLastSuccess(int lastSuccess) {
            FAILED.lastSuccess = lastSuccess;
            SUCCEEDED.lastSuccess = lastSuccess;
            AFTER_PERTURBATION.lastSuccess = lastSuccess;
        }
        public int lastSuccess() {
            return lastSuccess;
        }
    }

    public static enum mcmMode {RELAXATION, OPTIMIZATION}

    public MCM(TotalEnergy energy, ArrayList<Score> scoreFunctions, Score optimizationScore, Minimizer minimizer, Perturbation perturbation,
               TemperatureGenerator temperatureGenerator, int maxSteps) {
        super(energy, maxSteps, 1);
        this.energy = energy;
        this.minimizer = minimizer;
        this.temperatureGenerator = temperatureGenerator;
        this.perturbation = perturbation;
        this.scoreFunctions = scoreFunctions;
        this.optimizationScore = optimizationScore;
        for (AbstractEnergy e : energy.energyTerms()) {
            if (e instanceof TetherEnergy) {
                if (((TetherEnergy) e).allAtomsTether()) allAtomTether = (TetherEnergy) e;
            }
        }
        if (allAtomTether == null)
            throw new RuntimeException("Please add an all atoms Tether term for numerical stability at the end of the perturbation");
    }


    public OptimizerStatus run(Logger log) throws OptimizerException, UpdateableException, EvaluationException,AlignmentException {
        this.log = log;
        return run();
    }

    public OptimizerStatus run() throws OptimizerException, UpdateableException, EvaluationException,AlignmentException {
        double oldScore;
        double[][] oldCoordinates;
        double currentScore;
        double dS;

        Random randomNumberGenerator = MeshiProgram.randomNumberGenerator();

        energy.setCurrent();
        energy.on();
        allAtomTether.off();

        oldScore = -99999999; // first step always accepted
        System.out.println(" Old score" + oldScore);
        for (int step = 1; step <= maxSteps; step++) {
            energy.setCurrent();
            energy.on();
            ExcludedVol ev = (ExcludedVol)energy.getEnergyTerm(new ExcludedVol());
            ev.off();
            allAtomTether.off();
            double temperature = temperatureGenerator.next();
            energy.evaluate();
            System.out.println("\nxxxxxxxxxxxxxxxxx step # " + step + "; temperature " + temperature + " \ncurrent energy = \n" + energy.report(step));
            Utils.println("\n step # " + step + "; temperature " + temperature + " \ncurrent energy = \n" + energy.report(step));

            oldCoordinates = getOldCoordinates(energy);
            try {
                if (step > 1) perturb(perturbation, step);

                try {
                    energy.setCurrent();
                } catch (Exception ex) {
                    System.out.println(" This step failed.\ncould not setResidue the minimization energy to be the current one due to " + ex);
                    setCoordinates(energy, oldCoordinates);
                    continue;
                }
                energy.on();
                try {
                log.mcm(scoreFunctions, energy, step,mcmStepResult.AFTER_PERTURBATION);
                } catch(Exception ex){throw new RuntimeException("????????????????????");}
                // The great moment
                allAtomTether.reset();
                try {
                    minimizer.run();
                } catch (Exception ex) {
                    System.out.println(" This step failed.\nFirst (tethered) minimization failed due to " + ex);
                    ex.printStackTrace();
                    setCoordinates(energy, oldCoordinates);
                    continue;
                }
                allAtomTether.off();
                energy.evaluate();
                OptimizerStatus os = minimizer.run();

                System.out.println("MCM step " + os);
                /*if (log != null) {
                    try {
                        log.mcm(energy, step,mcmStepResult.FAILED);
                    }
                    catch (IOException ex) {
                        throw new RuntimeException(ex);
                    }
                }     */
                energy.on();
                allAtomTether.off();

//                for (Score scoreFunction:scoreFunctions)     {
//                    ((CombinedEnergyScore)scoreFunction).setEnergy(energy);
//                    scoreFunction.score(energy.energyInfo());
//                }
                currentScore = optimizationScore.score(energy.energyInfo());
                System.out.println("xxxxxxxxxxxxxxxxx oldScore = " + oldScore + "\n" + "currentScore = " + currentScore + "\n" + energy.report(999999));
                Utils.println("oldScore = " + oldScore + "\n" + "currentScore = " + currentScore + "\n" + energy.report(999999));
                dS = currentScore - oldScore;
                if (dS < 0) {
                    double rnd = randomNumberGenerator.nextDouble();
                    if (rnd > Math.exp(dS / temperature)) {//That is Metropolis criterion failed
                         if (log != null) {
                            try {
                                log.mcm(scoreFunctions, energy, step,mcmStepResult.FAILED);
                            } catch (IOException ex) {
                                throw new RuntimeException(ex);
                            }
                        }
                        setCoordinates(energy, oldCoordinates);
                        Utils.println("This step failed: dS = " + dS +
                                "  dS/temperature = " + dS / temperature +
                                "  Math.exp(dS/temperature) = " + Math.exp(dS / temperature) +
                                "  rnd = " + rnd);

                    } else {
                         if (log != null) {
                            try {
                                mcmStepResult.setLastSuccess(step);
                                log.mcm(scoreFunctions, energy, step,mcmStepResult.SUCCEEDED);
                                energy.distanceMatrix().writeNonBondedList(step+"_SUCCEEDED_1");
                            } catch (IOException ex) {
                                throw new RuntimeException(ex);
                            }
                        }
                        Utils.println("This step ded: dS = " + dS +
                                "  dS/temperature = " + dS / temperature +
                                "  Math.exp(dS/temperature) = " + Math.exp(dS / temperature) +
                                "  rnd = " + rnd);
                        oldScore = currentScore;
                    }
                } else {
                    System.out.println("This step succeeded: dS = " + dS);
                     if (log != null) {
                            try {
                                mcmStepResult.setLastSuccess(step);
                                energy.distanceMatrix().writeNonBondedList(step+"_SUCCEEDED_2");
                                log.mcm(scoreFunctions, energy, step,mcmStepResult.SUCCEEDED);
                            } catch (IOException ex) {
                                throw new RuntimeException(ex);
                            }
                     }
                    oldScore = currentScore;
                }
            }
            catch (OptimizerException ox) {
                energy.test();
                System.out.println("This step failed due to OptimizerException");
                setCoordinates(energy, oldCoordinates);
            }


        }
        energy.on();

        return OptimizerStatus.DONE;
    }

    public static double[][] getOldCoordinates(TotalEnergy energy) {
        double[][] coordinates = energy.coordinates();
        int length = coordinates.length;
        double[][] out = new double[length][2];

        for (int i = 0; i < length; i++) {
            out[i][0] = coordinates[i][0];
            out[i][1] = coordinates[i][1];
        }

        return out;
    }

    public static void setCoordinates(TotalEnergy energy, double[][] toSet) {
        double[][] coordinates = energy.coordinates();
        int length = coordinates.length;
        if (length != toSet.length) throw new RuntimeException("Weird parameters to MCM.setCoordinates");

        for (int i = 0; i < length; i++) {
            coordinates[i][0] = toSet[i][0];
            coordinates[i][1] = toSet[i][1];
        }
    }

    private static Perturbation[] toArray(Perturbation p) {
        Perturbation[] out = {p};
        return out;
    }

    protected void perturb(Perturbation perturbation, int step) throws OptimizerException, UpdateableException, EvaluationException{
         Minimizer.terminator.reset();
         perturbation.reset();
         AtomList atomListCopy = Utils.duplicateInAnewMolecularSystem(energy.atomList(),new MolecularSystem());
                    int counter = 0;
                    try {
                        perturbation.perturb();
                    }
                    catch (OptimizerException ox) {
                        System.out.println("A problem in perturbation #  " + counter);
                        System.out.println(ox);
                        ox.printStackTrace();
                        System.out.println("RMS reached\t" + counter + "\t" + step + "\t" + energy.atomList().getRms(atomListCopy));
                        if (ox.energy() != null) ox.energy().test();
                        System.out.println("\nContinuing\n");
                    }
                    System.out.println("RMS perturbation\t" + step + "\t" + energy.atomList().getRms(atomListCopy) + "\t" + log.rms());
        Minimizer.terminator.reset();
       }

}

/*
        InflateBySegments inflate = (InflateBySegments) energy.getEnergyTerm(new InflateBySegments());
        if (inflate == null) throw new RuntimeException("No point in MCM without inflate");
        if (iteration != 1) {
        inflate.isOn();
        System.out.println(optimizer.run());
        }
        inflate.off();
        System.out.println(optimizer.run());
*/
