/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

import meshi.energy.*;
import meshi.molecularElements.Protein;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.info.*;
/**

 */
public class ModelAnalyzer {
    public final Protein nativeStructure, model, originalModel;
    private TotalEnergy energy;
    private ResidueAlignmentMethod residueAlignmentMethod;

    public void setEnergy(TotalEnergy energy) {
        this.energy = energy;
    }

    public ModelAnalyzer(Protein model,
                         Protein nativeStructure,
                         Protein originalModel,
                         TotalEnergy energy,
                         ResidueAlignmentMethod residueAlignmentMethod) {
        this.nativeStructure = nativeStructure;
        this.model = model;
        this.originalModel = originalModel;
        this.energy = energy;
        this.residueAlignmentMethod = residueAlignmentMethod;
    }

    public double rms() throws Exception{
        return Rms.rms(nativeStructure, model, residueAlignmentMethod);

    }
    public double change() throws Exception{
        return Rms.rms(originalModel, model, residueAlignmentMethod);

    }

    public double rmsHeavy() throws Exception{
        return Rms.rmsHeavy(nativeStructure, model, residueAlignmentMethod);

    }

    public ProteinInfo analyze(String comment) throws UpdateableException, EvaluationException,AlignmentException {
        ProteinInfo out = new ProteinInfo(comment, model.sourceFile().getAbsolutePath(), model.name(), model);
        if (nativeStructure != null) {
            double[] gdt_ts = Rms.gdt(nativeStructure, model);
            out.add(new DoubleInfoElement(InfoType.GDT_TS, "GDT_TS form native structure", gdt_ts[0]));
            double[] originalGdt = Rms.gdt(nativeStructure, originalModel);
            out.add(new DoubleInfoElement(InfoType.DELTA_GDT_TS, "delta GDT_TS with respect to the original unrefined", gdt_ts[0] - originalGdt[0]));
            double[] gdt_ha = Rms.gdt(nativeStructure, model,GDTcalculator.Type.HA);
            out.add(new DoubleInfoElement(InfoType.GDT_HA, "GDT_HA form native structure", gdt_ha[0]));
            try {
                double[] originalGdtHa = Rms.gdt(nativeStructure, originalModel,GDTcalculator.Type.HA);
                out.add(new DoubleInfoElement(InfoType.DELTA_GDT_HA, "delta GDT_HA with respect to the original unrefined", gdt_ha[0] - originalGdtHa[0]));
                double rms = Rms.rms(nativeStructure, model, residueAlignmentMethod);
                out.add(new DoubleInfoElement(InfoType.RMS, "RMS form native structure", rms));
                double originalRms = Rms.rms(nativeStructure, originalModel, residueAlignmentMethod);
                out.add(new DoubleInfoElement(InfoType.DELTA_RMS, "delta RMS with respect to the original unrefined model", rms - originalRms));
                double rmsHeavy = Rms.rmsHeavy(nativeStructure, model, residueAlignmentMethod);
                out.add(new DoubleInfoElement(InfoType.RMS_HEAVY, "Heavy atoms RMS form native structure", rmsHeavy));
                double originalRmsHeavy = Rms.rmsHeavy(nativeStructure, originalModel, residueAlignmentMethod);
                out.add(new DoubleInfoElement(InfoType.DELTA_RMS_HEAVY, "delta RMS of heavy atoms with respect to the original unrefined model", rmsHeavy - originalRmsHeavy));
            }catch (Exception ex ) {throw new RuntimeException(ex.getMessage());}
            double rmsFromOriginal;
            try {
                rmsFromOriginal = Rms.rms(originalModel, model, residueAlignmentMethod);
        out.add(new DoubleInfoElement(InfoType.CHANGE, "structural change (in RMS) due to refinement", rmsFromOriginal));
            } catch (Exception ex ) {throw new RuntimeException(ex.getMessage());}
        }

        if (energy != null) {
            double e = energy.evaluateAll();
            out.add(new DoubleInfoElement(InfoType.ENERGY, "energy", e));
            for (AbstractEnergy energyElement : energy.energyTerms()) {
                out.add(new EnergyInfoElement(energyElement.info()));
                for (DoubleInfoElement element : energyElement.info().infoList()) {
                    out.add(new DoubleInfoElement(element));
                }
            }
        }
        return out;
    }
}

