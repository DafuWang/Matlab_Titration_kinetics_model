function [bulk,SO4_Ba_Cl]=ion_z_C0_G_lambda_infinity_bulk(ion)

bulk.z=[ion.K.z, ion.Na.z, ion.Ca.z, ion.H.z, ion.OH.z, ion.NO3.z, ion.Cl.z];
bulk.C0=[ion.K.C0_bulk, ion.Na.C0_bulk, ion.Ca.C0_bulk, ion.H.C0_bulk, ion.OH.C0_bulk, ion.NO3.C0_bulk, ion.Cl.C0_bulk];
bulk.G=[ion.K.G, ion.Na.G, ion.Ca.G, ion.H.G, ion.OH.G, ion.NO3.G, ion.Cl.G];
bulk.lambda_infinity=[ion.K.lambda_infinity, ion.Na.lambda_infinity, ion.Ca.lambda_infinity, ion.H.lambda_infinity, ion.OH.lambda_infinity, ion.NO3.lambda_infinity, ion.Cl.lambda_infinity];

SO4_Ba_Cl.z=[ion.SO4.z, ion.Ba.z, ion.Cl.z];
SO4_Ba_Cl.C0=[ion.SO4.C0_Na2SO4,ion.Ba.C0_BaCl2,ion.Cl.C0_BaCl2];
SO4_Ba_Cl.G=[ion.SO4.G, ion.Ba.G, ion.Cl.G];
SO4_Ba_Cl.lambda_infinity=[ion.SO4.lambda_infinity, ion.Ba.lambda_infinity, ion.Cl.lambda_infinity];

