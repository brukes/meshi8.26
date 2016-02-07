/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.scoringFunctions;

import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.UpdateableException;
import meshi.util.Utils;
import meshi.util.file.MeshiLineReader;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.DoubleInfoElementList;
import meshi.util.info.InfoType;
import org.xml.sax.SAXException;

import javax.xml.parsers.ParserConfigurationException;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 10/11/2010
 * Time: 23:18:00
 * To change this template use File | Settings | File Templates.
 */
public class CombinedEnergyScore implements Score, KeyWords {
    private NormalizedEnergyScore[] scoreFunctions;
    private double interdecileRange = -9999;
//    private TotalEnergy energy;
    private double[] sumWeights;
    private DoubleInfoElement lengthElement, changeElement;
    private String name;
    private String[] fields;

    public CombinedEnergyScore(String parametersFileName, String name) throws IOException, UpdateableException,EvaluationException,
                                                                                         ParserConfigurationException,SAXException{
//        this.energy = null;
        this.name = name;
        ConfigurationArray configurationArray= new ConfigurationArray(new File(parametersFileName));
        scoreFunctions = new NormalizedEnergyScore[configurationArray.size()];
        int iConfig = 0;
        fields = configurationArray.get(0).fields;
        for (Configuration configuration:configurationArray)     {
            scoreFunctions[iConfig] = new NormalizedEnergyScore(configuration);
            iConfig++;
        }
        sumWeights = new double[configurationArray.size()];
        lengthElement = null;
        changeElement = null;
    }

    public String toString() {
        return name;
    }

//    public void setEnergy(TotalEnergy energy) {
//        this.energy = energy;
//    }
    public void setLengthElement(int length){
        lengthElement = new DoubleInfoElement(InfoType.LENGTH,"The length of the protein",length);
    }
    public void setChangeElement(double change) {
        changeElement = new DoubleInfoElement(InfoType.CHANGE, "structural change (in RMS) due to refinement", change);
    }

    public boolean fieldRequired(String s) {
        for (String field :fields)
            if (field.equals(s)) return true;
        return false;
    }
    public double score(DoubleInfoElementList energyInfo) throws UpdateableException, EvaluationException{
//        if (energy == null)
//            throw new RuntimeException("energy is null");

        if (fieldRequired("length")) {
            if (!energyInfo.contains(lengthElement)) {
                if (lengthElement == null)
                    throw new RuntimeException("score function "+toString()+": lengthElement is null");
                energyInfo.addElement(lengthElement);
            }

        }

        if (fieldRequired("change")) {
            if (!energyInfo.contains(changeElement)) {
                if (changeElement == null)
                    throw new RuntimeException("changeElement is null");
                energyInfo.addElement(changeElement);
            }
        }



        for (int i = scoreFunctions.length-1; i >= 0; i--) {// Makes it easier to compare to the MATLAB version
            NormalizedEnergyScore scoreFunction = scoreFunctions[i];
            if (scoreFunction == null)
                throw new RuntimeException("This is weird "+i);
            scoreFunction.calcScore(energyInfo);
        }
        if (scoreFunctions.length == 1) {
            interdecileRange = 0;
            return scoreFunctions[0].getScore();
        }
        Arrays.sort(scoreFunctions);
        double sum = 0;
        for (int i = 0; i < sumWeights.length; i++)
            sumWeights[i] = sum = sum + scoreFunctions[i].getWeight();
        int midPoint = Arrays.binarySearch(sumWeights,sum/2);
        if (midPoint < 0 ) //sum/2 is not explicitly in the array
            midPoint = -midPoint;
        if (sumWeights.length >= 10) {
            int low  = Arrays.binarySearch(sumWeights,sum/10);
            if (low < 0) low = -low;
            int high = Arrays.binarySearch(sumWeights,9*sum/10);
            if (high < 0) high = -high;
            interdecileRange = scoreFunctions[high].getScore()-scoreFunctions[low].getScore();
        }
        else interdecileRange = 0;

        int topMedianIndex = midPoint; //It is the index above  sum/2
        double topMedian = scoreFunctions[topMedianIndex].getScore();
        int bottomMedianIndex = midPoint-1;
        double bottomWeight = sumWeights[topMedianIndex]-sum/2;
        double topWeight;
        double bottomMedian;
        if (topMedianIndex == 0) {
            bottomMedian = 0;
            topWeight = sumWeights[topMedianIndex];
        }
        else {
            bottomMedian = scoreFunctions[bottomMedianIndex].getScore();
            topWeight = sum/2-sumWeights[bottomMedianIndex];
        }

        return (topMedian*topWeight+bottomMedian*bottomWeight)/(topWeight+bottomWeight);

        //return (scoreFunctions[midPoint].getScore()+scoreFunctions[midPoint-1].getScore())/2;


    }

    public double interdecileRange() {
        double out = interdecileRange;
        interdecileRange = -9999;
        return out;
    }
}

