/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.pairwiseNonBondedTerms.CooperativeSumma;

import meshi.energy.EnergyInfoElement;
import meshi.util.info.InfoType;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 17/11/2010
 * Time: 00:30:13
 * To change this template use File | Settings | File Templates.
 */
public class CooperativeSummaInfoElement extends EnergyInfoElement {
    public final EnergyInfoElement COOPERATIVE_SUMMA_POLAR;
    public final EnergyInfoElement COOPERATIVE_SUMMA_NON_POLAR;
    public final EnergyInfoElement COOPERATIVE_SUMMA_NEUTRAL;
    public final EnergyInfoElement COOPERATIVE_SUMMA_POLAR_NN_OO;
    public final EnergyInfoElement COOPERATIVE_SUMMA_POLAR_BB;

    public CooperativeSummaInfoElement(InfoType type, String comment, EnergyInfoElement polar,
                                       EnergyInfoElement nonPolar,
                                       EnergyInfoElement neutral,
                                       EnergyInfoElement polarNNOO,
                                       EnergyInfoElement polarBB,
                                       double weight) {
        super(type, comment, weight);
        infoList().add(COOPERATIVE_SUMMA_POLAR = polar);
        infoList().add(COOPERATIVE_SUMMA_NON_POLAR = nonPolar);
        infoList().add(COOPERATIVE_SUMMA_NEUTRAL = neutral);
        infoList().add(COOPERATIVE_SUMMA_POLAR_NN_OO = polarNNOO);
        infoList().add(COOPERATIVE_SUMMA_POLAR_BB = polarBB);
     }
}
