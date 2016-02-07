/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.scoringFunctions;

import meshi.energy.TotalEnergy;
import meshi.molecularElements.Protein;
import meshi.util.KeyWords;
import meshi.util.info.*;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 10/11/2010
 * Time: 23:18:00
 * To change this template use File | Settings | File Templates.
 */
public class NormalizedEnergyScore implements Comparable {
    private Sigmoid sigmoid;
    private Configuration configuration;
    private double score, weight;



    public NormalizedEnergyScore(Configuration configuration)  {
        sigmoid = new Sigmoid(configuration);
        this.configuration = configuration;

    }

    public void calcScore(DoubleInfoElementList energyInfo)  {
        double[] fieldValues      = getValues(energyInfo, configuration.fields);
        double[] weightedValues   = new double[configuration.fields.length];
        double[] normalizerValues = getValues(energyInfo, configuration.normalizers);
        double tempScore = 0;
        for (int iField = 0; iField < fieldValues.length; iField++) {
            double normalizationFactor = getNormalizationFactor(normalizerValues,configuration.exponentIndices[iField]);
            weightedValues[iField] = fieldValues[iField] *
                                     configuration.coefs[iField] *
                                     normalizationFactor;
            tempScore += weightedValues[iField];
            //System.out.println("yyyyyyy "+configuration.fields[iField]+" "+fieldValues[iField]+" "+configuration.coefs[iField]+" "+normalizationFactor+" "+tempScore);
        }
//        double tempScore = 0;
//        for (double d:weightedValues)
//            tempScore += d;
        score = sigmoid.backSigmoid(tempScore);
        weight = sigmoid.derivative(score);
    }

    public double getScore() {
        return score;
    }
    public double getWeight() {
        return weight;
    }


    private double getNormalizationFactor(double[] normalizersValues, int[] exponentsIndices){
        double factor = 1;
        for (int iNormalizer = 0; iNormalizer < normalizersValues.length; iNormalizer++){
            try {
                factor *= normalizersExponent(normalizersValues[iNormalizer],exponentsIndices[iNormalizer]);
            }
            catch (RuntimeException ex) {
                System.out.println("error in getNormalizationFactor\niNormalizer = "+iNormalizer);
                throw ex;
            }
        }
        return factor;
    }

    private double normalizersExponent(double normalizerValue, int exponentIndex){
        // Note that the indices here are taken from a MATLAB program and they (MATLAB people) never heard of index 0.
        if (normalizerValue < 0)
            throw new RuntimeException("Base cannot be smaller than 0, but it is "+normalizerValue);
        if (normalizerValue == 0) exponentIndex = 2;
        if (exponentIndex-1 >= configuration.exponentValues.length) return 0;
        if (configuration.exponentValues[exponentIndex-1] == 1)    return normalizerValue;
        if (configuration.exponentValues[exponentIndex-1] == 0.5)  return Math.sqrt(normalizerValue);
        if (configuration.exponentValues[exponentIndex-1] == 0)    return 1;
        if (configuration.exponentValues[exponentIndex-1] == -0.5) return 1/Math.sqrt(normalizerValue);
        if (configuration.exponentValues[exponentIndex-1] == -1)   return 1/normalizerValue;
        if (configuration.exponentValues[exponentIndex-1] == -2)   return 1/(normalizerValue*normalizerValue);
        if (configuration.exponentValues[exponentIndex-1] == -3)   return 1/(normalizerValue*normalizerValue*normalizerValue);
        throw new RuntimeException("weird value "+configuration.exponentValues[exponentIndex]);

    }

        private double[] getValues(DoubleInfoElementList energyInfo, String[] fields){
            double[] out = new double[fields.length];
            for (int iField = 0; iField < fields.length; iField++) {
                DoubleInfoElement element = null;
                for (DoubleInfoElement currentElement:energyInfo)
                    if (currentElement.type.tag.equals(fields[iField])){
                        if (element != null)
                            throw new RuntimeException("field "+fields[iField]+" appears more than once "+element+" "+currentElement);
                        else
                            element = currentElement;
                    }
                if (element == null)  {
                    for (DoubleInfoElement currentElement:energyInfo)
                        System.out.println(currentElement.type.tag);
                    throw new MissingFieldException("field "+fields[iField]+" is missing.");
                }
                out[iField] = element.value();
            }

            return out;
        }

    public int compareTo(Object obj) {
        if (score - ((NormalizedEnergyScore) obj).score >0) return 1;
        if (score - ((NormalizedEnergyScore) obj).score <0) return -1;
        return 0;
    }
}

/*
        double[] normalizersValues = new double configuration;
        if (energyInfo.size() != coefs.length)
                                throw new RuntimeException("This is weird energyInfo.size() = "+energyInfo.size()+"  weights.length = "+coefs.length);
        double out = 0;
//        System.out.println("yyyyyyy "+normalizers.size());
        for (int iInfoType = 0; iInfoType < energyInfo.size(); iInfoType++) {
            if (coefs[iInfoType] != 0) {
//                System.out.println("xxxxxxxxx "+weights[iInfoType]+" "+lengthExponents[iInfoType]+" "+normalizersExponents[iInfoType].get(0).doubleValue()+" "+energyInfo.get(iInfoType));
                double s =   energyInfo.get(iInfoType).doubleValue();
                if (s >= VERY_LARGE) {
                   throw new RuntimeException("Score error: "+energyInfo.get(iInfoType)+" has a "+s+"  value, which is probably an error. Maybe it failed to find its parameters file.") ;
                }
                s *= coefs[iInfoType];
                if (lengthExponents[iInfoType] <-99)
                    s = 0;
                else {
                    s *= myPower(protein.chain().numberOfNonDummyResidues(), lengthExponents[iInfoType]);
                    for (int iNormalizer = 0; iNormalizer<normalizers.size(); iNormalizer++) {
                        double   exponent   = normalizersExponents[iInfoType].get(iNormalizer).doubleValue();
                        if (exponent < -99) s = 0;
                        else s *= myPower(normalizersValues[iNormalizer],exponent);
                    }
                }
                out += s;
            }

        }
        return out;
    }

    private double myPower(double base,double exponent) {
        if (base <= 0) throw new RuntimeException("Base cannot be smaller than 0");
        if (exponent == 1)    return base;
        if (exponent == 0.5)  return Math.sqrt(base);
        if (exponent == 0)    return 1;
        if (exponent == -0.5) return 1/Math.sqrt(base);
        if (exponent == -1)   return 1/base;
        if (exponent == -2)   return 1/(base*base);
        if (exponent == -3)   return 1/(base*base*base);
        throw new RuntimeException("weird exponent "+exponent);
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

    public double[] getNormalizersValues(DoubleInfoElementList energyInfo) {
        double[] out = new double[normalizers.size()];

        for (int i = 0; i<normalizers.size(); i++){
            InfoType infoType = normalizers.get(i);
            boolean found = false;
            for (DoubleInfoElement element : energyInfo) {
                  if (element.type == infoType) out[i] = element.doubleValue();
                  found = true;
            }
            if (!found) throw new RuntimeException("This is weird. Failed to find "+infoType);
        }

        return out;
    }

  public static MeshiInfoElementList getWeights(ArrayList<String> lines, ArrayList<InfoType> normalizers) {
    MeshiInfoElementList out;
    out = new MeshiInfoElementList();
    for (String line : lines) {
        ArrayList<String> words = new ArrayList<String>();
        StringTokenizer tokenizer = new StringTokenizer(line);
        while(tokenizer.hasMoreTokens())
            words.add(tokenizer.nextToken());
        if (words.size() != normalizers.size()+3 )
            throw new RuntimeException("This is weird");
        InfoType infoType = InfoType.getTypeByTag(words.get(0));
        if (infoType != null) {
            ScoreInfoElement scoreInfoElement = new ScoreInfoElement(infoType,
                                         "Score weight of " + infoType,
                                         Double.parseDouble(words.get(1)),
                                         Double.parseDouble(words.get(2)));
            for (int i = 3; i < words.size();i++)
                scoreInfoElement.addNormalizationExponent(words.get(i));
            out.add(scoreInfoElement);
        }
    }
    return out;
  }

    private static class ScoreInfoElement extends DoubleInfoElement {
        public double lengthExponent;
        public ArrayList<Double> normalizersExponents = new ArrayList<Double>();
        public ScoreInfoElement(InfoType type, String comment, double weight, double lengthExponent) {
            super(type, comment, weight);
            this.lengthExponent = lengthExponent;
        }
        public void addNormalizationExponent(String word) {
                normalizersExponents.add(new Double(word));
        }
    }


}
  */