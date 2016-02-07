/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.pairwiseNonBondedTerms.CooperativeSumma;

import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.AtomCore;
import meshi.parameters.AtomType;


/**
 * Created by IntelliJ IDEA.
 * User: TM
 * Date: Feb 15, 2010
 * Time: 10:46:36 PM
 * To change this template use File | Settings | File Templates.
 */
public class CooperativeZSummaPolarElement_simplified {
    private static final double Z_ALPHA = 10, Z_WEIGHT = 0.0000000001;
    protected final static int NumOfClusters = CooperativeZ5typesParameters.NumOfClusters;
    protected final static double[] expOfZScoredAtomEnergies = new double[NumOfClusters];
    protected final static double[] zScoredAtomEnergiesDerivatives = new double[NumOfClusters];
    protected final static double[] zScoredAtomEnergies = new double[NumOfClusters];
    protected final static double[] zPenalty = new double[NumOfClusters];
    protected final static double[] zPenaltyDerivative = new double[NumOfClusters];
    private double[] zEnergy, zEnergyDerivative;
    private double allButI;
    private int sizeOfMS;
    private final double z_weight;

    private static double zMinusZ_ALPHA;
    double summPlusEpsilon, summ, summDivSummPlusEpsilon;
    double[][] mean = new double[NumOfClusters][AtomType.values().length];
    double[][] std = new double[NumOfClusters][AtomType.values().length];
    double[][] proportion = new double[NumOfClusters][AtomType.values().length];
    double[] atomEnergies;
    int nAtoms;
    int initialPickNumber = 0;  //starts from 0 to 2
    double koef = 1;
    double EPSILON = 0.00001;

    static final double WIDTH = CooperativeSumma5PolarTypesEnergy.WIDTH;
    static final double HEIHT = CooperativeSumma5PolarTypesEnergy.HEIHT;
    static double M = (HEIHT + 1) * WIDTH;
    double weight;

    double[][] derivedExpOfZScoredAtomEnergies;
    double[] z = new double[NumOfClusters];
    double[] energyForEachND = new double[NumOfClusters],
            derivativeForEachND = new double[NumOfClusters];
    double energy;
    String name;

    CooperativeZSummaPolarElement_simplified(double[][] mean, double[][] std, double[][] proportion, int nAtoms, double[] atomEnergies, double allWeight, DistanceMatrix distanceMatrix, String name, double myWeight) {
        this.mean = mean;
        this.std = std;
        this.proportion = proportion;
        this.nAtoms = nAtoms;
        sizeOfMS = distanceMatrix.molecularSystem.size();
        if (nAtoms == 0)
            throw new RuntimeException("Zero number of atoms in " + name + this + ". Weird polarElement definition.");
        this.atomEnergies = atomEnergies;
        this.weight = allWeight * myWeight;
        z_weight = Z_WEIGHT * myWeight;
        this.name = name;
        derivedExpOfZScoredAtomEnergies = new double[NumOfClusters][sizeOfMS];
        zEnergy = new double[sizeOfMS];
        zEnergyDerivative = new double[sizeOfMS];


    }

    void resetZ() {
        for (int i = initialPickNumber; i < NumOfClusters; i++) {
            z[i] = 0;
            energyForEachND[i] = 0;
            derivativeForEachND[i] = 0;
            for (int k = 0; k < derivedExpOfZScoredAtomEnergies[i].length; k++)
                derivedExpOfZScoredAtomEnergies[i][k] = 0;
        }
        energy = 0;
    }

    void addZScore(AtomCore atom) {
        // Comments by Chen
        // 1. zScoredAtomEnergies atomEnergies are transformed to Zscore - a vector of three numbers for the three peaks.
        //    The highest values for the furthest peak.
        // 2. expOfZScoredAtomEnergies[i] - yet another transformation on the three numbers - values between 1 (just below a peak) to zero (far away).
        // 3. summ - the sum of the three expOfZScoredAtomEnergies[i]
        summ = 0;
        int atomNumber = atom.number;
        int type = atom.type().ordinal();
        for (int i = initialPickNumber; i < NumOfClusters; i++) {
            zScoredAtomEnergies[i] = zScore(atomEnergies[atomNumber], mean[i][type], std[i][type]);
            expOfZScoredAtomEnergies[i] = proportion[i][type] * Math.exp(-0.5 * koef * (zScoredAtomEnergies[i] * zScoredAtomEnergies[i]));
            summ += expOfZScoredAtomEnergies[i];
        }

        // Comments by Chen
        // expEnergy - normalizer
        // zArray [i] - sum over all atomOne of the contribution of the ith Gaussian
        double expEnergy;
        for (int i = initialPickNumber; i < NumOfClusters; i++) {
            expEnergy = expOfZScoredAtomEnergies[i] / summPlusEpsilon;
            z[i] += (zScoredAtomEnergies[i] * expEnergy);
              }

    }


    double zScore(double curE, double mean, double std) {
        if (std == 0)
            //               return 0;
            throw new RuntimeException("Zero std value in " + this);
        return (curE - mean) / std;
    }

    void makeEnergy() {
        for (int i = initialPickNumber; i < NumOfClusters; i++) {
            z[i] /= nAtoms;
            calcEnergy2(i);
            energy += energyForEachND[i];
        }

    }

    void calcEnergy2(int i) {
        /* ------zArray^2 */
        double zz = z[i];
        double z2 = zz * zz;
        energyForEachND[i] = weight * (M / WIDTH - (M + z2) / (z2 + WIDTH));
    }


    void evaluateAtoms(AtomCore atom1, AtomCore atom2, SummaAttribute summaAttribute, int order) {
        double derivative, force;
        for (int i = initialPickNumber; i < NumOfClusters; i++) {
            derivative = derivedExpOfZScoredAtomEnergies[i][atom2.number];
            //       if (derivative == 0.0)
            //            // throw new RuntimeException("Something is weird in Parameters "+derivative+name);
            //    System.out.println("Something is weird in Parameters "+derivative+name);

            force = derivativeForEachND[i] * derivative / 2;
            force = force * order;
            if (!atom1.status().frozen())
                atom1.addForce(force * summaAttribute.fx, force * summaAttribute.fy, force * summaAttribute.fz);
            if (!atom2.status().frozen())
                atom2.addForce(-force * summaAttribute.fx, -force * summaAttribute.fy, -force * summaAttribute.fz);
        }
        force = zEnergyDerivative[atom2.number] * order / 2.0;

        if ((!(force < 0)) & (!(force == 0)) & (!(force > 0)))
            throw new RuntimeException("this is weird1 " + force + " " + atom2);
        if (!atom1.status().frozen())
            atom1.addForce(force * summaAttribute.fx, force * summaAttribute.fy, force * summaAttribute.fz);
        if (!atom2.status().frozen())
            atom2.addForce(-force * summaAttribute.fx, -force * summaAttribute.fy, -force * summaAttribute.fz);

    }

}
