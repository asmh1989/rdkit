use std::fmt::{Debug, Formatter};

use cxx::let_cxx_string;
use rdkit_sys::*;

use crate::{Fingerprint, RWMol};

pub struct ROMol {
    pub(crate) ptr: cxx::SharedPtr<ro_mol_ffi::ROMol>,
}

impl ROMol {
    pub fn from_smile(smile: &str) -> Result<Self, cxx::Exception> {
        let_cxx_string!(smile_cxx_string = smile);
        let ptr = ro_mol_ffi::smiles_to_mol(&smile_cxx_string)?;
        Ok(Self { ptr })
    }

    pub fn from_smile_with_params(
        smile: &str,
        params: &SmilesParserParams,
    ) -> Result<Self, cxx::Exception> {
        let_cxx_string!(smile_cxx_string = smile);
        let ptr = ro_mol_ffi::smiles_to_mol_with_params(&smile_cxx_string, params.ptr.clone())?;
        Ok(Self { ptr })
    }

    pub fn as_smile(&self) -> String {
        ro_mol_ffi::mol_to_smiles(self.ptr.clone())
    }

    pub fn as_rw_mol(&self, quick_copy: bool, conf_id: i32) -> RWMol {
        let ptr = rdkit_sys::rw_mol_ffi::rw_mol_from_ro_mol(self.ptr.clone(), quick_copy, conf_id);
        RWMol { ptr }
    }

    pub fn remove_hs(&self) -> ROMol {
        let ptr = rdkit_sys::ro_mol_ffi::remove_hs(self.ptr.clone());
        ROMol { ptr }
    }

    pub fn fingerprint(&self) -> Fingerprint {
        let ptr = fingerprint_ffi::fingerprint_mol(self.ptr.clone());
        Fingerprint::new(ptr)
    }

    pub fn fingerprint_2_vec(&self) -> Vec<String> {
        let fingerprint = fingerprint_ffi::fingerprint_mol2(self.ptr.clone());
        let bytes: Vec<String> = fingerprint.into_iter().map(|x| (*x).to_string()).collect();
        bytes
    }
}

impl Debug for ROMol {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let smile = self.as_smile();
        f.debug_tuple("ROMol").field(&smile).finish()
    }
}

impl Clone for ROMol {
    fn clone(&self) -> Self {
        ROMol {
            ptr: rdkit_sys::ro_mol_ffi::copy_mol(self.ptr.clone()),
        }
    }
}

pub struct SmilesParserParams {
    pub(crate) ptr: cxx::SharedPtr<ro_mol_ffi::SmilesParserParams>,
}

impl SmilesParserParams {
    pub fn sanitize(&mut self, value: bool) {
        rdkit_sys::ro_mol_ffi::smiles_parser_params_set_sanitize(self.ptr.clone(), value);
    }
}

impl Default for SmilesParserParams {
    fn default() -> Self {
        SmilesParserParams {
            ptr: rdkit_sys::ro_mol_ffi::new_smiles_parser_params(),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use crate::Properties;

    use super::*;

    #[test]
    fn test_name() {
        let mol = ROMol::from_smile("c1ccccc1C(=O)NC").unwrap();
        let properties = Properties::new();
        let computed: HashMap<String, f64> = properties.compute_properties(&mol);
        assert_eq!(*computed.get("NumAtoms").unwrap(), 19.0);
    }
}
