/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.scoringFunctions;

import meshi.energy.EvaluationException;
import meshi.util.UpdateableException;
import meshi.util.info.DoubleInfoElementList;

/**

 */
public interface Score  {
    public double score(DoubleInfoElementList energyInfo) throws UpdateableException, EvaluationException;


}
