/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.pairwiseNonBondedTerms;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.geometry.Distance;
import meshi.geometry.DistanceLists;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.Atom;
import meshi.util.Updateable;

import java.util.ArrayList;

public abstract class NonBondedEnergyTerm extends AbstractEnergy {
    protected DistanceMatrix distanceMatrix;
    /**
     * Must be initialized isOn subclass construction.
     */
    protected NonBondedEnergyElement energyElement = null;

    public NonBondedEnergyTerm() {
        super();
    }

    /**
     * Creates a new <code>NonBondedEnergyTerm</code> instance.
     *
     * @param updateableResources an <code>Object[]</code> value
     * @param distanceMatrix      a <code>DistanceMatrix</code> value
     */
    public NonBondedEnergyTerm(Updateable[] updateableResources,
                               EnergyInfoElement info,
                               DistanceMatrix distanceMatrix) {
        super(updateableResources, info);
        this.distanceMatrix = distanceMatrix;
    }


    /**
     * Evaluates energy for each distance
     *
     * @return a sum of all energy elements
     */
    public EnergyInfoElement evaluate()throws EvaluationException {

        double energy = 0;
        DistanceLists nonBondedList = distanceMatrix.nonBondedList();
        for (ArrayList<Distance> row :nonBondedList) {
            for (Distance distance : row) {
                if (!distance.mode().frozen) {
                    energyElement.set(distance);
                    energy += energyElement.evaluate();
                }
            }
        }
        info.setValue(energy);
        return info;
    }

    /**
     * Describe <code>evaluateAtoms</code> method here.
     */
    public void evaluateAtoms()throws EvaluationException {
        if (on) {
            DistanceLists nonBondedList = distanceMatrix.nonBondedList();
            for (ArrayList<Distance> row :nonBondedList) {
                for (Distance distance : row) {
                    if (!distance.mode().frozen) {
                        energyElement.set(distance);
                        energyElement.evaluateAtoms();
                    }
                }
            }
        }
    }

    /**
     * Testing of one atomOne in all energy elements
     *
     * @param totalEnergy a <code>TotalEnergy</code> value
     * @param atom        an criminal <code>Atom</code> value
     */
    public void test(TotalEnergy totalEnergy, Atom atom) throws EvaluationException{
        if (energyElement == null) throw new RuntimeException("energyElement is null");
        if (!on) System.out.println("" + this + " is off");
        DistanceLists nonBondedList = distanceMatrix.nonBondedList();
        for (ArrayList<Distance> row :nonBondedList) {
            for (Distance distance : row) {
                if (distance == null) {
                    //System.out.println("Null distance in NonBondedList");
                    System.out.println("Weird nonBondedList:");
                    for (Distance d : row) System.out.println(d);
                    throw new RuntimeException("Null distance in NonBondedList");
                }
                if ((distance.atom1() == null) || (distance.atom2() == null)) {
                    //System.out.println("Null distance in NonBondedList"+nonBonded);
                    System.out.println("Weird nonBondedList:");
                    for (Distance d : row) System.out.println(d);
                    throw new RuntimeException("Null atomOne in distance " + distance);
                }

                if (distance == null) throw new RuntimeException("nonBonded is null");
                energyElement.set(distance);
                energyElement.test(totalEnergy, atom);
            }
        }
    }
}
