package meshi.energy.secondaryStructure;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EnergyInfoElement;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.info.InfoType;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 20/05/12
 * Time: 11:34
 * To change this template use File | Settings | File Templates.
 */
public class SecondaryStructureFractionCreator extends EnergyCreator{
    public SecondaryStructureFractionCreator() {
        super(InfoType.SECONDARY_STRUCTURE_FRACTION);
    }
    EnergyInfoElement info = new EnergyInfoElement(InfoType.SECONDARY_STRUCTURE_FRACTION,"The fraction of secondary structure residues.");
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
                                                    CommandList commands) {
        term = new SecondaryStructureFraction(info,protein);
        return term;
    }
}
