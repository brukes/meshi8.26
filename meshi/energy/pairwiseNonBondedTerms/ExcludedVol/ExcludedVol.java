/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

/**
 *Excluded Volume potential.
 *
 *A potential that creates a strong repulsion between atoms when they are getting
 *near their VDW radius. There is no attraction part like in the VDW.
 *The functional form of the term is:  
 *
 *dis =                EV
 *
 *[0,sigma]  C*(dis-sigma)^4
 *[sigma,Inf]          0
 *
 *width is setResidue in ExcludedVolParameters.java. Currently it is 0.2 Ang.
 *width is the transition zone (in Ang) where the energy change in the forth power 0.0 to 1.0.
 *
 **/
package meshi.energy.pairwiseNonBondedTerms.ExcludedVol;

import meshi.energy.EnergyInfoElement;
import meshi.energy.EvaluationException;
import meshi.energy.pairwiseNonBondedTerms.NonBondedEnergyTerm;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.DistanceLists;
import meshi.geometry.Distance;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomCore;
import meshi.util.Classes;
import meshi.util.Utils;
import meshi.util.filters.Filter;

public class ExcludedVol extends NonBondedEnergyTerm implements Classes {
    private final Filter filter;
    private final int excludedVolumType;
    private final double rFactor;
    private final double rMax;
    private final ExcludedVolParametersList parametersList;

    public ExcludedVol(){filter = null; excludedVolumType = -1; rFactor = -1; rMax = -1; parametersList = null;}
    public ExcludedVol(DistanceMatrix distanceMatrix,
                       ExcludedVolParametersList parametersList,
                       int type,
                       EnergyInfoElement info,
                       double rfac,
                       Filter filter) {
        super(toArray(distanceMatrix), info, distanceMatrix);
        comment = "ExcludedVol";
        energyElement = new ExcludedVolEnergyElement(parametersList, distanceMatrix, type, weight, rfac);
        this.filter = filter;
        this.rFactor = rfac;
        this.excludedVolumType = type;
        this.parametersList = parametersList;
        rMax = distanceMatrix.rMax();
    }


    public void scaleWeight(double factor) {weight = weight*factor;}
    public EnergyInfoElement evaluate() throws EvaluationException{
        ExcludedVolParameters parameters;
        Atom atom1,atom2;
        double dEdD;
        double DMinusSig;
        double D3;
        double localSigma;
        double sigma,C1,C2;
        int atom1Type, atom2Type; //well, their ordinals
        int atom1ResidueNumber;
        double dDx,dDy,dDz;
        double atom1Dx, atom1Dy,atom1Dz;
        boolean atom1IsBackbone;
        if (!on) {
            info.setValue(0);
            return info;
        }
        double energy = 0;
        DistanceLists nonBondedList;
        if (filter == null)
            nonBondedList = distanceMatrix.nonBondedList();
        else {
            nonBondedList = new DistanceLists(100);
            Utils.filter(distanceMatrix.nonBondedList(), filter, nonBondedList);
        }
        if (nonBondedList == null)
            System.out.println("!!!!!!");

        for (DistanceList distanceRow : nonBondedList) {
            atom1 = distanceRow.atomOne.atom;
            atom1Type = atom1.type().ordinal();
            atom1ResidueNumber = atom1.residueNumber();
            atom1IsBackbone = atom1.isBackbone();
            atom1Dx = atom1Dy = atom1Dz = 0;
            for (Distance distance : distanceRow) {
                dDx = distance.dDistanceDx();
                dDy = distance.dDistanceDy();
                dDz = distance.dDistanceDz();

                atom2 = distance.atom2();
                    atom2Type = atom2.type().ordinal();
                    if (atom1Type > atom2Type)
                        parameters = parametersList.parameters(atom1Type,atom2Type);
                    else
                        parameters = parametersList.parameters(atom2Type,atom1Type);
                    sigma = parameters.sigma;
                    C1 = parameters.C1;
                    C2 = parameters.C2;
                    if (sigma > rMax)
                        throw new RuntimeException("Excluded Volume: sigma=" + sigma + " and it is larger " +
                            "than rMax=" + rMax);
                    localSigma = sigma;
                if (rFactor < 0.99) {
                    if (atom1ResidueNumber == atom2.residueNumber())
                        localSigma = 0;
                    else
                        localSigma = sigma * rFactor;
                }
                double dis = distance.distance();
                if (dis > localSigma) continue;

                DMinusSig = (dis - localSigma);
                switch (excludedVolumType) {
                    case 0:
                        if (atom1IsBackbone || atom2.isBackbone()) {
                            D3 = C2 * DMinusSig * DMinusSig * DMinusSig * weight;
                            energy += D3 * DMinusSig;
                            dEdD = 4 * D3;
                        } else {
                            D3 = C1 * DMinusSig * weight;
                            energy += D3 * DMinusSig;
                            dEdD = 2 * D3;
                        }
                        break;
                    case 1:
                        if (atom1IsBackbone && atom2.isBackbone()) {
                            D3 = 1.0 / (0.1 * 0.1) * DMinusSig * weight;
                            energy += D3 * DMinusSig;
                            dEdD = 2 * D3;
                        } else {
                            energy += 0.0;
                            dEdD = 0.0;
                        }
                        break;
                    case 2:
                        if (atom1IsBackbone && atom2.isBackbone()) {
                            D3 = DMinusSig * weight;
                            energy += D3 * DMinusSig;
                            dEdD = 2 * D3;
                        } else {
                            energy += 0.0;
                            dEdD = 0.0;
                        }
                        break;
                    default:
                        throw new RuntimeException("An unknown excluded volume type");
                }
                if (!atom1.frozen())  {
                    atom1Dx += -1 * dEdD * dDx;
                    atom1Dy += -1 * dEdD * dDy;
                    atom1Dz += -1 * dEdD * dDz;
                }
                if (!atom2.frozen())
                    atom2.core.addForce(dEdD * dDx,
                                        dEdD * dDy,
                                        dEdD * dDz);
            }

            atom1.core.addForce(atom1Dx, atom1Dy, atom1Dz);
        }
        info.setValue(energy);
        return info;
    }//evaluate
}

	
 
