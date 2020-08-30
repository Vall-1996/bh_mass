# bh_mass
The 'critical scripts' are:

- ozdes_spectra.py (which relies on ozdes_tools.py) mostly for just looking at the spectra.

- Janie's getPhoto.py (https://github.com/jhoormann/OzDES_getPhoto) for fetching magnitudes from photometry,
  OR querying the databases yourself for the magnitudes (which is what I ended up doing).
  - for DES DR1 queries, I used the DES easyaccess package (https://github.com/mgckind/easyaccess) and readQuery.py to read the output. 
  - for NOAO queries, I used their online query interface (https://datalab.noao.edu/query.php) and pulled from the image headers. 
  
- Janie's calibSpec.py (https://github.com/jhoormann/OzDES_calibSpec) for calibrating the spectra using the magnitudes from the photometry.

- fit_spectra.py (which relies on generic_peak.py) for fitting Gaussians to emission lines. This is an MCMC fitter and fits up to three peaks at once.

- mass.py (which relies on ozdes_mass_tools.py) for calculating BH masses from fit parameters.
  - H beta and H alpha constants come from this paper: https://arxiv.org/pdf/astro-ph/0508335.pdf
  - Mg II constants come from John's paper: https://arxiv.org/pdf/1410.1539.pdf

Everything else is 'old scripts' or 'temp scripts' that I decided to toss in here as a back up.

- ozdes_legacy.py, ozdes_database.py: old / scrapped scripts from summer 2019 when I started this project.
- FWHM.py, integrate.py, straighten.py: temp scripts when I was figuring out how to calculate mass from fit parameters.
- read_obs.py, all_mass.py: temp scripts I made for file-reading convenience and cranking through calculations.
- spec-xxxx-xxxxx-xxxx.fits: some SDSS spectra I practiced fitting on, from this paper: https://arxiv.org/pdf/1509.03634.pdf
- links.txt: just a backup of some links I had bookmarked as I had to refer back to them a lot.

Stuff backed up on irulan:

- The OzDES spectra as SpARCS_xxxx.fits and SVA1_COADD_xxxx.fits, under the folder 'spectra'.
- The DES DR1 queries as DESDR_xxxx.tab, under the folder 'desdr1_queries'.
- The NOAO queries are NOT backed up, as I deleted each image after scraping the header for magnitudes (space limitations).
- The DES DR1 magnitudes as mag_dr1.txt, under the folder 'magnitudes'.
- The NOAO-queried magnitudes as mag_noao.txt, under the folder 'magnitudes'.
- The OzDES F star spectra as OzDES-DR2_xxxxx.fits and FStars.fits, under the folder 'f_stars'.
- A copy of this repository, under the folder 'scripts'.
- The emission line fit parameters and black hole masses as bh_mass_results.xlsx, under the folder 'fits_and_masses'.


