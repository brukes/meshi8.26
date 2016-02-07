/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.scoringFunctions;

import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.molecularElements.Protein;
import meshi.util.*;
import meshi.util.info.*;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 10/11/2010
 * Time: 23:18:00
 * To change this template use File | Settings | File Templates.
 */
public class EnergyScore implements Score, KeyWords {
    private double[] weights;
    private TotalEnergy energy;
    public enum Normalization {NOT_NORMALIZED, NORMALIZED_BY_LENGTH, IGNORED};
    private Normalization[] normalizations;
    private Protein protein;
    private MeshiInfoElementList energyScoreWeights;
    public static final double VERY_LARGE = Double.MAX_VALUE/2;



    public EnergyScore(TotalEnergy energy,
                       MeshiInfoElementList energyScoreWeights,
                       Protein protein) throws UpdateableException,EvaluationException {
        this.energy = energy;
        this.energyScoreWeights = energyScoreWeights;
        if (energyScoreWeights != null) {
            DoubleInfoElementList energyInfo = energy.energyInfo();
            weights = new double[energyInfo.size()];
            normalizations = new Normalization[energyInfo.size()];
            for (int i = 0; i < weights.length; i++) weights[i] = 0;
            for (MeshiInfoElement weightElement : energyScoreWeights) {
                int index = getIndex(energyInfo, weightElement);
                if (index != -1) {
                    //System.out.println("Found energyInfo for "+weightElement); // DEBUG
                    weights[index] = weightElement.doubleValue();
                    normalizations[index] = ((ScoreInfoElement) weightElement).normalization;
                }
                else {
                    throw new RuntimeException("Cannot find "+weightElement+" in "+energyInfo);
                }
            }
            for (int i = 0; i < normalizations.length; i++) {
                Normalization normalization = normalizations[i];
                if (normalization == null){
                    weights[i] = 0;
                    normalizations[i] =Normalization.IGNORED;
                }
            }
        }
        this.protein = protein;
    }

    public double score(DoubleInfoElementList energyInfo) throws UpdateableException,EvaluationException {
        if (energyScoreWeights == null) {
              return -1*energy.evaluate();
        }
        if (energyInfo.size() != weights.length)
                                throw new RuntimeException("This is weird energyInfo.size() = "+energyInfo.size()+"  weights.length = "+weights.length);
        double out = 0;
        for (int i = 0; i < energyInfo.size(); i++) {
            if (weights[i] != 0) {
                double s =   energyInfo.get(i).doubleValue();
                if (s >= VERY_LARGE) {
                    //throw new RuntimeException("Score error: "+energyInfo.get(i)+" has a "+s+"  value, which is probably an error. Maybe it failed to find its parameters file.") ;
                   s = 0;// Probably  missing file due to lack of homologs.
                }
                s *= weights[i];
                if (normalizations[i].equals(Normalization.NORMALIZED_BY_LENGTH))
                    s = s/ protein.chain().numberOfNonDummyResidues();
                out += s;
            }
        }
        return out;
    }

    private static int getIndex(MeshiInfoElementList list, MeshiInfoElement weightElement) {
        for (int i = 0; i < list.size(); i++) {
            System.out.println("getIndex2: Looking for type "+weightElement.type); // DEBUG
            if (list.get(i).type == weightElement.type) {
                if (weightElement.type.valueType != InfoValueType.DOUBLE)
                    throw new RuntimeException(" Weird energy weight element " + weightElement);
                return i;
            }
        }
        return -1;
    }

    private static int getIndex(DoubleInfoElementList list, MeshiInfoElement weightElement) {
        for (int i = 0; i < list.size(); i++) {
            //System.out.println("getIndex1: Looking for type "+weightElement.type+" vs "+list.get(i).type); // DEBUG
            if (list.get(i).type == weightElement.type) {
                //System.out.println("      getIndex1: Found type "+list.get(i).type); 
                if (weightElement.type.valueType != InfoValueType.DOUBLE)
                    throw new RuntimeException(" Weird energy weight element " + weightElement);
                return i;
            }
        }
        return -1;
    }

    public static MeshiInfoElementList getWeights(CommandList commandList) {
        return getWeights(commandList, SCORE_WEIGHT);
    }
        public static MeshiInfoElementList getWeights(CommandList commandList, Key key) {
    MeshiInfoElementList out = null;
        if(commandList.keyExists(key)) {
                out = new MeshiInfoElementList();
                CommandList list = commandList.firstWordFilter(key);
                for (Command command : list) {
                        InfoType infoType = InfoType.getTypeByTag(command.secondWord());
                        if (infoType != null) {
                                out.add(new ScoreInfoElement(infoType, "Score weight of " + infoType, command.thirdWordDouble(),command.fourthWord()));
                                //System.out.println("Score weight of " + infoType+" tag="+command.secondWord());
                        }
                }
        }
        return out;
    }

    private static class ScoreInfoElement extends DoubleInfoElement {
        public Normalization normalization = null;
        public ScoreInfoElement(InfoType type, String comment, double value, String flag) {
            super(type, comment, value);
            for (Normalization norm : Normalization.values()) {
                if (norm.toString().equals(flag)) normalization = norm;
            }
            if (normalization == null)
                throw new RuntimeException("Cannot find normalization mode for "+this);
        }
    }

}
