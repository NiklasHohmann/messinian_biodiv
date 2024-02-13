# messinian_biodiv

Code for "Late Miocene transformation of Mediterranean Sea biodiversity"  
Project webpage: [REMARE project](https://sites.google.com/view/kagiadi/projects/remare)

## Authors

__Niklas Hohmann__  (creator and maintainer of repository)
Utrecht University  
email: n.h.hohmann [at] uu.nl  
Web page: [www.uu.nl/staff/NHohmann](https://www.uu.nl/staff/NHHohmann)  
ORCID: [0000-0003-1559-1838](https://orcid.org/0000-0003-1559-1838)

__Konstantina Agiadi__ (principal investigator)  
University of Vienna  
email: konstantina.agiadi [at] univie.ac.at  
Web page: [sites.google.com/view/kagiadi](https://sites.google.com/view/kagiadi)  
ORCID: [0000-0001-8073-559X](https://orcid.org/0000-0001-8073-559X)  

## License

Apache 2.0, see LICENSE file for full text.

## Requirements

Base R (version >= 4) and the RStudio IDE.

## Reproducing Results

In the RStudio IDE, open the file _messinian_biodiv.Rproj_. This opens the RProject of the same name. Then, run

```R
source("code/analysis.R")
```

to reproduce the results. This will (1) produce all figures in the _figs_ folder and (2) generate variables `sr_change` and `eco_ind_median` in your workspace that contain median values of species richness and ecological indices. To inspect them, use

```R
sr_change
eco_ind_median
```

## Repository Structure

* _code_ : folder with code
  * _analysis.R_ : code for main analysis
  * _helper_function.R_ : aux functions (select data from DB, rarefy species richness & other ecological indices)
* _data_ : folder containing the Messinina Database
* _figs_ : folder for figures. Initially all subfolders are empty, they will be filled after the code is run (see section _Reproducing Results_)
  * _eco_timeslice_comp_ : figs of comparison of ecol indices
  * _eco_timeslice_comp_regional_ : figs of comparison of ecol indices per region
  * _sr_through_time_ : figs of species richness through time
  * _sr_through_time_regional_ : figs of species richness through time per region
* _.gitignore_ : untracked files
* _LICENSE_ : Apache 2.0 license text
* _messinian_biodiv.RProj_ : Rproject file
* _README_ : README file

## Funding

This work was supported by the Austrian Science Fund (FWF) project “Late Miocene Mediterranean Marine Ecosystem Crisis” (2022–2026), Grant DOI 10.55776/V986 (PI: K.Agiadi)."  
