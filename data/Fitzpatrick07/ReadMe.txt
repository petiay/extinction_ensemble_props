J/ApJ/663/320       IR-through-UV extinction curve          (Fitzpatrick+, 2007)
================================================================================
An analysis of the shapes of interstellar extinction curves.
V. The IR-through-UV curve morphology.
    Fitzpatrick E.L., Massa D.
   <Astrophys. J., 663, 320-341 (2007)>
   =2007ApJ...663..320F
================================================================================
ADC_Keywords: Models ; Extinction
Keywords: dust, extinction - methods: data analysis

Abstract:
    We study the IR-through-UV wavelength dependence of 328 Galactic
    interstellar extinction curves affecting normal, near-main-sequence B
    and late O stars. We derive the curves using a new technique that
    employs stellar atmosphere models in lieu of unreddened "standard"
    stars. Under ideal conditions, this technique is capable of virtually
    eliminating spectral mismatch errors in the curves. In general, it
    lends itself to a quantitative assessment of the errors and enables a
    rigorous testing of the significance of relationships between various
    curve parameters, regardless of whether their uncertainties are
    correlated.

File Summary:
--------------------------------------------------------------------------------
 FileName   Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe         80        .   This file
table1.dat     67      328   Basic Data for Survey Stars
table3.dat     87      328   Best-Fit Parameters for Survey Stars
table4.dat    138      328   Best-Fit Extinction Curve Parameters for
                              Survey Stars
refs.dat       64       83   References
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table1.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1- 22  A22   ---     Name      Star name (1)
  24- 38  A15   ---     SpType    MK Spectral type (2)
  40- 44  F5.2  mag     Vmag      The V band magnitude
  46- 49  I4    pc      Dist      ? Heliocentric distance (3)
  51- 57  F7.3  deg     GLON      Galactic longitude
  59- 64  F6.2  deg     GLAT      Galactic latitude
  66- 67  I2    ---     Ref       ? Reference, in refs.dat file
--------------------------------------------------------------------------------

Note (1): The stars are listed in order of increasing Right Ascension using
     the most commonly adopted forms of their names. The first preference
     was ``HDnnn'', followed by ``BDnnn'', etc. There are 185 survey stars
     which are members of open clusters or associations. The identity of
     the cluster or association is either contained in the star name itself
     (e.g., NGC 457 Pesch 34) or is given in parentheses after the star's
     name.

Note (2): Spectral types were selected from those given in the SIMBAD
     database, and the source of the adopted types is shown in the "Ref".
     When multiple types were available for a particular star, we selected
     one based on our own preferred ranking of the sources.

Note (3): Distances are:
    NGC 2244 = distance is from Perez et al. (1987PASP...99.1050P);
    NGC 3293 = distance is from Balona & Crampton (1974MNRAS.166..203B);
    Trumpler 14 and 16 = distances are from Massey & Johnson
          (1993, Cat. J/AJ/105/980);
    Cep OB3 = distance is from Crawford & Barnes (1970AJ.....75..952C).

     The distances to all other clusters or associations are from the Open
     Clusters and Galactic Structure database maintained by Wilton S. Dias,
     Jacques Lepine, Bruno S. Alessi, and Andre Moitinho, Cat. B/ocl, and
     http://www.astro.iag.usp.br/~wilton/ .

     For the non-cluster stars, distances were calculated using the E(B-V)
     values from this study and the absolute magnitudes from Turner
     (1980ApJ...240..137T) (for mid-B and earlier types) and Blaauw 1963
     (in Basic Astronomical Data, ed. K. A. Strand (Chicago:
     Univ. Chicago Press), chap. 20) (for mid-B and later types).
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table3.dat
--------------------------------------------------------------------------------
   Bytes Format Units     Label     Explanations
--------------------------------------------------------------------------------
   1- 22  A22   ---       Name      Star name
  24- 28  I5    K         Teff      Model effective temperature (2)
  30- 33  I4    K       e_Teff      1{sigma} uncertainty in Teff
  35- 38  F4.2  [cm/s2]   log(g)    Log of model surface gravity (3)
  40- 43  F4.2  [cm/s2] e_log(g)    1{sigma} uncertainty in log(g)
  45- 49  F5.2  [Sun]     [m/H]     Log of model metallicity (4)
  51- 54  F4.2  [Sun]   e_[m/H]     ? 1{sigma} uncertainty in [m/H]
  56- 59  F4.1  km/s      Vturb     Model turbulent velocity (5)
  61- 63  F3.1  km/s    e_Vturb     ? 1{sigma} uncertainty in Vturb
  65- 70  F6.4  mas       theta     Model angular radius
  72- 77  F6.4  mas     e_theta     1{sigma} uncertainty in theta
  79- 82  F4.2  mag       E(B-V)    Model reddening
  84- 87  F4.2  mag     e_E(B-V)    1{sigma} uncertainty in E(B-V)
--------------------------------------------------------------------------------
Note (2): For the O stars analyzed using the TLUSTY atmosphere models, the
     values of Teff were adopted from the Spectral Type vs. T_eff_ relation
     given in Table 2. These stars can be identified by their 1-{sigma}
     uncertainties, which are +/-1000K.

     Table 2: Adopted Temperature Scale for Main-Sequence O Stars
     -------------------
     SpType    Teff (K)
     -------------------
     O6        40000
     O6.5      38500
     O7        37000
     O7.5      36500
     O8        36000
     O8.5      34750
     O9        33500
     O9.5      32750
     B0        32000
     -------------------

Note (3): For stars in clusters, the surface gravities are determined as
     discussed in Sect. 3.1 and rely on stellar evolution models and
     cluster distance determinations. Surface gravities for non-cluster
     stars are not always well-determined, because of a lack of specific
     spectroscopic indicators. In some cases, the best-fit solutions for
     these stars indicated physically unlikely results (i.e., log(g)>~4.3
     or log(g)<~3.0). For these stars, a value of log(g)=3.9 was assumed
     (which is the mean log(g) of the rest of the sample) and a 1-{sigma}
     uncertainty of +/-0.2 was incorporated in the error analysis. These
     cases can be identified by log(g) entries of "3.9+/-0.2".

Note (4): For the O stars in the sample, our fitting procedure utilized
     solar abundance TLUSTY models. For these stars the values of [m/H] 
     are indicated by entries of "0" without uncertainties.

Note (5): For the O stars, the adopted TLUSTY models incorporate
     Vturb=10km/s. For these stars the values of Vturb are indicated by
     entries of "10" without uncertainties.

     For the B stars, which were modeled using ATLAS9 models, the values of
     Vturb were determined by the fitting procedure, but were constrained
     to lie between 0 and 10km/s.

     Stars whose best-fit SED models required these limiting values are
     indicated by Vturb entries of "0" or "10", without error bars. The
     uncertainties for stars with best-fit Vturb values close to these
     limits may be underestimated due to this truncation.
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table4.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1- 22  A22   ---     Name      Star name
  24- 28  F5.3  ---     x0        The UV x_0_ coefficient (2)
  30- 34  F5.3  ---   e_x0        1-{sigma} uncertainty in x0
  36- 39  F4.2  ---     gamma     UV {gamma} coefficient (2)
  41- 44  F4.2  ---   e_gamma     1-{sigma} uncertainty in {gamma}
  46- 50  F5.2  ---     c1        UV c_1_ coefficient (2)
  52- 55  F4.2  ---   e_c1        ? 1-{sigma} uncertainty in c1
  57- 61  F5.2  ---     c2        UV c_2_ coefficient (2)
  63- 66  F4.2  ---   e_c2        1-{sigma} uncertainty in c2
  68- 72  F5.2  ---     c3        UV c_3_ coefficient (2)
  74- 77  F4.2  ---   e_c3        1-{sigma} uncertainty in c3
  79- 82  F4.2  ---     c4        UV c_4_ coefficient (2)
  84- 87  F4.2  ---   e_c4        1-{sigma} uncertainty in c4
  89- 92  F4.2  ---     c5        UV c_5_ coefficient (2)
  94- 97  F4.2  ---   e_c5        1-{sigma} uncertainty in c5
  99-102  F4.2  ---     O1        ? Optical spline O_1_ point (3)
 104-107  F4.2  ---   e_O1        ? 1-{sigma} uncertainty in O1
 109-112  F4.2  ---     O2        Optical spline O_2_ point (3)
 114-118  F5.2  ---     O3        Optical spline O_3_ point (3)
 120-123  F4.2  ---     R(V)      IR R(V) coefficient (4)
 125-128  F4.2  ---   e_R(V)      ? 1-{sigma} uncertainty in R(V)
 130-133  F4.2  ---     kIR       IR k_IR_ coefficient (4)
 135-138  F4.2  ---   e_kIR       ? 1-{sigma} uncertainty in kIR
--------------------------------------------------------------------------------
Note (2): Extinction curve is defined by
       k({lambda}-V) = c1+c2x+c3.D(x,x0,{gamma}) for x<=c5 and
                     = c1+c2x+c3.D(x,x0,{gamma}) = c4(x-c5)^2^ for x>c5.
     where x={lambda}^-1^, in units of inverse microns (um-1) and
     D(x,x0,{gamma})=x^2^/[(x^2^-x0^2^)+x^2^{gamma}^2^]

     For the stars HD237019, HD18352, and HD25443 the long wavelength IUE
     spectra are incomplete. For these cases we constrained the UV linear
     extinction component to follow the relation c1=2.18-2.91*c2 from
     Fitzpatrick (2004, in ASP Conf. Ser. 309, 33). For these stars we 
     list uncertainties for the c2 values but not for the c1 values.

Note (3): The uncertainties in the O2 and O3 optical spline points (at
     wavelengths of 4000 and 5530 Angstroms, respectively) are typically
     0.01 or less and are not listed. For several stars, those without U
     band photometry, we did not solve for the O1 point at 3300 Angstroms.

Note (4): R(V) is the ratio of reddening to extinction at V. For field
     stars without IR photometry, we assumed R(V)=3.1 and kIR=1.11, with
     the latter based on the relation kIR=0.63*R(V)-0.84 from Fitzpatrick
     (2004, in ASP Conf. Ser. 309, 33) clusters, we adopted the mean R(V)
     of the other cluster members and a value of kIR based on the
     aforementioned relation. These assumed values are listed in the Table
     without uncertainties. Several survey stars have apparently noisy JHK
     data and yielded very uncertain values of kIR. For these, we
     ultimately derived the extinction curve by solving for the best-fit
     value of R(V) with kIR constrained to follow the Fitzpatrick (2004, in
     ASP Conf. Ser. 309, 33) relation. The resultant R(V) values are listed
     with their uncertainties while the kIR values are listed without
     uncertainties.
--------------------------------------------------------------------------------

Byte-by-byte Description of file: refs.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1-  2  I2    ---     Ref       Reference number
   4- 22  A19   ---     BibCode   BibCode
  24- 45  A22   ---     Aut       Author's name
  46- 64  A19   ---     Com       Comments
--------------------------------------------------------------------------------

History:
    From electronic version of the journal

References:
   Fitzpatrick & Massa,  Paper I     1986ApJ...307..286F
   Fitzpatrick & Massa,  Paper II    1988ApJ...328..734F
   Fitzpatrick & Massa,  Paper III   1990ApJS...72..163F
   Fitzpatrick & Massa,  Paper IV    2005AJ....130.1127F
   Fitzpatrick & Massa,  Paper VI    2009ApJ...699.1209F
================================================================================
(End)                  Greg Schwarz [AAS], Patricia Vannier [CDS]    14-Aug-2009
