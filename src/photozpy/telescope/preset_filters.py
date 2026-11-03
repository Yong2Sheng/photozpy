from .filterset import FilterSet

SDSS_PRIMED = FilterSet(
    filterset_name="SDSS_PRIMED",
    standard_to_alias={
        "SDSS_u'": ("Sloan u'2", "SDSS u", "Sloan_u'2", "SDSS_u", "SDSS_up"),
        "SDSS_g'": ("Sloan g'2", "SDSS g", "Sloan_g'2", "SDSS_g", "SDSS_gp"),
        "SDSS_r'": ("Sloan r'2", "SDSS r", "Sloan_r'2", "SDSS_r", "SDSS_rp"),
        "SDSS_i'": ("Sloan i'2", "SDSS i", "Sloan_i'2", "SDSS_i", "SDSS_ip"),
        "SDSS_z'": ("Sloan z'2", "SDSS z", "Sloan_z'2", "SDSS_z", "SDSS_zp"),
    },
)


BESSEL = FilterSet(
    filterset_name="BESSEL",
    standard_to_alias={
        "Bessel_U": ("U",),
        "Bessel_B": ("B",),
        "Bessel_V": ("V",),
        "Bessel_I": ("I",),
        "Bessel_R": ("R",),
    },
)
