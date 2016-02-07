/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.rg;

import meshi.energy.EnergyInfoElement;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.InfoType;

public class RgInfoElement extends EnergyInfoElement {


    public final double all, RGlogaritmicWeight, RGratioWeight, RGpolarWeight, RGnonPolarWeight, RGbackboneWeight;
    public final DoubleInfoElement N_RG_HSS, E_HSS_HCOIL, E_HSS_BSS, E_HSS_BCOIL, E_HSS_PSS, E_HSS_PCOIL, HSS, BSS, PSS, HSS_HCOIL, HSS_BSS, HSS_BCOIL, HSS_PSS, HSS_PCOIL;

    public RgInfoElement(double all, double RGlogaritmicWeight, double RGratioWeight, double RGpolarWeight, double RGnonPolarWeight, double RGbackboneWeight) {
        super(InfoType.RG, "Radius of Gyration related terms ");
        this.all = all;
        this.RGlogaritmicWeight = RGlogaritmicWeight;
        this.RGratioWeight = RGratioWeight;
        this.RGpolarWeight = RGpolarWeight;
        this.RGnonPolarWeight = RGnonPolarWeight;
        this.RGbackboneWeight = RGbackboneWeight;
        N_RG_HSS = infoList().addElement(new DoubleInfoElement(InfoType.N_RG_HSS, "log(Length)  Vs. log(RG of hydrophobic side-chains in secondary structure elements", Double.MIN_VALUE));
        E_HSS_HCOIL = infoList().addElement(new DoubleInfoElement(InfoType.E_HSS_HCOIL, "log(RG of hydrophobic side-chains in secondary structure elements Vs of hydrophobic side-chains in coil", Double.MIN_VALUE));
        E_HSS_BSS = infoList().addElement(new DoubleInfoElement(InfoType.E_HSS_BSS, "log(RG of hydrophobic side-chains in secondary structure elements Vs CAs in SSEs", Double.MIN_VALUE));
        E_HSS_BCOIL = infoList().addElement(new DoubleInfoElement(InfoType.E_HSS_BCOIL, "log(RG of hydrophobic side-chains in secondary structure elements Vs CAs in coil", Double.MIN_VALUE));
        E_HSS_PSS = infoList().addElement(new DoubleInfoElement(InfoType.E_HSS_CSS, "log(RG of hydrophobic side-chains in secondary structure elements Vs Polar side chain atoms in SSEs", Double.MIN_VALUE));
        E_HSS_PCOIL = infoList().addElement(new DoubleInfoElement(InfoType.E_HSS_CCOIL, "log(RG of hydrophobic side-chains in secondary structure elements Vs Polar side chain atoms in Coil", Double.MIN_VALUE));
        HSS = infoList().addElement(new DoubleInfoElement(InfoType.HSS, "log(RG of hydrophobic side-chains in secondary structure elements", Double.MIN_VALUE));
        BSS = infoList().addElement(new DoubleInfoElement(InfoType.BSS, "log(RG of CAs in secondary structure elements", Double.MIN_VALUE));
        PSS = infoList().addElement(new DoubleInfoElement(InfoType.CSS, "log(RG of polar atoms in secondary structure elements", Double.MIN_VALUE));
        HSS_HCOIL = infoList().addElement(new DoubleInfoElement(InfoType.HSS_HCOIL, "log(RG of hydrophobic side-chains in secondary structure elements Vs of hydrophobic side-chains in coil", Double.MIN_VALUE));
        HSS_BSS = infoList().addElement(new DoubleInfoElement(InfoType.HSS_BSS, "log(RG of hydrophobic side-chains in secondary structure elements Vs CAs in SSEs", Double.MIN_VALUE));
        HSS_BCOIL = infoList().addElement(new DoubleInfoElement(InfoType.HSS_BCOIL, "log(RG of hydrophobic side-chains in secondary structure elements Vs CAs in coil", Double.MIN_VALUE));
        HSS_PSS = infoList().addElement(new DoubleInfoElement(InfoType.HSS_CSS, "log(RG of hydrophobic side-chains in secondary structure elements Vs Polar side chain atoms in SSEs", Double.MIN_VALUE));
        HSS_PCOIL = infoList().addElement(new DoubleInfoElement(InfoType.HSS_CCOIL, "log(RG of hydrophobic side-chains in secondary structure elements Vs Polar side chain atoms in Coil", Double.MIN_VALUE));

    }
}