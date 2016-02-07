/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvate;

import meshi.energy.EnergyInfoElement;
import meshi.util.info.InfoType;

/**

 */

public class SolvateInfoElement extends EnergyInfoElement {
    public final EnergyInfoElement scPolarInfo, bbPolarInfo, scCarbonInfo, stdInfo;

    public SolvateInfoElement(InfoType type, String comment, double weightSCPolarSolvate,
                              double weightSCCarbonSolvate,
                              double weightBBPolarSolvate) {
        super(type, comment, 1);
        scPolarInfo = (EnergyInfoElement) infoList().addElement(new EnergyInfoElement(InfoType.SOLVATION_SC_POLAR, "Solvation term for sidechain polar atoms", weightSCPolarSolvate));
        scCarbonInfo = (EnergyInfoElement) infoList().addElement(new EnergyInfoElement(InfoType.SOLVATION_SC_CARBON, "Solvation term for sidechain carbon atoms", weightSCCarbonSolvate));
        bbPolarInfo = (EnergyInfoElement) infoList().addElement(new EnergyInfoElement(InfoType.SOLVATION_BB_POLAR, "Solvation term for backbone  polar atoms", weightBBPolarSolvate));
        stdInfo     =  (EnergyInfoElement) infoList().addElement(new EnergyInfoElement(InfoType.SOLVATION_STD," Distribution of solvation energy over the protein"));
    }

    EnergyInfoElement stdElememt() {
        return stdInfo;
    }
}
