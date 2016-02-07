/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvate;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.info.InfoType;

/**
 * The different between this creator and the "SolvateCreatorRegularHB" (see the latter main comment for a lot of details), is in the
 * way the hydrogen-bond distance sigmoid is setResidue. Here parameters are chosen so that the hydrogen-bond starts its drop to 0 at a
 * donor-acceptor distance of 3.1 Ang, and is 0 only for donor-acceptor distance of 4.0 Ang. This is good to increase the basin of
 * attraction of direct minimizations protocols.
 */

public class SolvateCreatorHBforMinimization extends EnergyCreator {

    // The different weights relevent to this term.
    private double weightSCPolarSolvate = -9999.0;
    private double weightSCCarbonSolvate = -9999.0;
    private double weightBBPolarSolvate = -9999.0;
    private double weightHB = -9999.0;
    private static SolvateParametersList parametersList = null;
    private String comment = null;
    private DistanceMatrix distanceMatrix = null;

    public SolvateCreatorHBforMinimization(double weightSCPolarSolvate, double weightSCCarbonSolvate,
                                           double weightBBPolarSolvate, double weightHB) {
        this(weightSCPolarSolvate, weightSCCarbonSolvate, weightBBPolarSolvate, weightHB, "Solvation");
    }

    public SolvateCreatorHBforMinimization(double weightSCPolarSolvate, double weightSCCarbonSolvate,
                                           double weightBBPolarSolvate, double weightHB, String comment) {
        super(InfoType.SOLVATION);
        this.weightSCPolarSolvate = weightSCPolarSolvate;
        this.weightSCCarbonSolvate = weightSCCarbonSolvate;
        this.weightBBPolarSolvate = weightBBPolarSolvate;
        this.weightHB = weightHB;
        this.comment = comment;
    }

    public SolvateCreatorHBforMinimization(double weight) {
        super(weight);
        this.weightSCPolarSolvate = weight;
        this.weightSCCarbonSolvate = weight;
        this.weightBBPolarSolvate = weight;
        this.weightHB = weight;
    }

    public SolvateCreatorHBforMinimization() {
        super(InfoType.SOLVATION);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
                                           CommandList commands) {
        if ((term != null) && (this.distanceMatrix == distanceMatrix)) //chen 3/3/12
                return term;
            this.distanceMatrix = distanceMatrix;
            if (weightSCPolarSolvate < 0) {
            weightSCPolarSolvate = commands.firstWordFilter(InfoType.SOLVATION).secondWord(InfoType.SOLVATION_SC_POLAR.tag).thirdWordDouble();
            weightSCCarbonSolvate = commands.firstWordFilter(InfoType.SOLVATION).secondWord(InfoType.SOLVATION_SC_CARBON.tag).thirdWordDouble();
            weightBBPolarSolvate = commands.firstWordFilter(InfoType.SOLVATION).secondWord(InfoType.SOLVATION_BB_POLAR.tag).thirdWordDouble();
            weightHB = commands.firstWordFilter(InfoType.SOLVATION).secondWord(InfoType.SOLVATION_HB.tag).thirdWordDouble();
        }
        if (parametersList == null) {
            // Appanding the path to the list of parameter filenames.
            String[] strlist = new String[SOLVATE_MINIMIZE_HB_PARAMETERS.length];
            String pathname = parametersDirectory(commands).concat("/");
            for (int cc = 0; cc < SOLVATE_MINIMIZE_HB_PARAMETERS.length; cc++)
                strlist[cc] = pathname.concat(SOLVATE_MINIMIZE_HB_PARAMETERS[cc]);
            parametersList = new SolvateParametersList(strlist);
        }
        weightSCPolarSolvate *= weight();
        weightSCCarbonSolvate *= weight();
        weightBBPolarSolvate *= weight();
        weightHB *= weight();

        SolvateInfoElement info = new SolvateInfoElement(InfoType.SOLVATION, "Solvation energy",
                weightSCPolarSolvate,
                weightSCCarbonSolvate,
                weightBBPolarSolvate);

        return term = new SolvateEnergy(protein.atoms(),
                distanceMatrix,
                parametersList,
                info);
    }

}
