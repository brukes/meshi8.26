/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.tether;

import meshi.energy.simpleEnergyTerms.*;

import meshi.energy.EnergyCreator;

import java.util.*;

import meshi.util.*;
import meshi.util.filters.*;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.*;

import meshi.geometry.DistanceMatrix;

/**
 * Created by IntelliJ IDEA.
 * User: guykarl
 * Date: 23/06/2005
 * Time: 12:32:27
 * To change this template use File | Settings | File Templates.
 */

import meshi.energy.*;

import meshi.molecularElements.*;
import meshi.util.info.InfoType;


public class TetherCreator extends SimpleEnergyCreator implements KeyWords {
    String comment = null;

    public TetherCreator(InfoType type) {
        super(type);
    }


    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
                                           CommandList commands) {
        if (term != null) return term;
        AtomList atomList = protein.atoms();
        double[][] initialLocationsMatrix = new double[3][atomList.size()];
        for (Iterator atoms = protein.atoms().iterator(); atoms.hasNext();) {
            Atom atom = (Atom) atoms.next();
            if (!atom.nowhere()) {
                int number = atom.number();
                try {
                    initialLocationsMatrix[0][number] = atom.x();
                    initialLocationsMatrix[1][number] = atom.y();
                    initialLocationsMatrix[2][number] = atom.z();
                }
                catch (RuntimeException ex) {
                    System.out.println("Problem in TetherCreator. createEnergyTerm while processing atomOne number :" + number + "\n" + atom);
                    throw ex;
                }
            }
        }
        EnergyInfoElement info = new EnergyInfoElement(infoType, "Tether energy ", weight);
        term = new TetherEnergy(atomList, initialLocationsMatrix, info, comment);
        return term;
    }
}

