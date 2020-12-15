"""
Test case for the Kuehn et a (2020) GA subduction model.
Test tables generated using R function as published at
github page
"""
from openquake.hazardlib.gsim.kuehn_2020 import (
    KuehnEtAl2020SInter,
    KuehnEtAl2020SInterAlaska,
    KuehnEtAl2020SInterCascadia,
    KuehnEtAl2020SInterCentralAmericaMexico,
    KuehnEtAl2020SInterJapan,
    KuehnEtAl2020SInterNewZealand,
    KuehnEtAl2020SInterSouthAmerica,
    KuehnEtAl2020SInterTaiwan,
    KuehnEtAl2020SSlab,
    KuehnEtAl2020SSlabAlaska,
    KuehnEtAl2020SSlabCascadia,
    KuehnEtAl2020SSlabCentralAmericaMexico,
    KuehnEtAl2020SSlabJapan,
    KuehnEtAl2020SSlabNewZealand,
    KuehnEtAl2020SSlabSouthAmerica,
    KuehnEtAl2020SSlabTaiwan)
from openquake.hazardlib.tests.gsim.utils import BaseGSIMTestCase

# Interface
class KuehnEtAl2020SInterTestCase(BaseGSIMTestCase):
    GSIM_CLASS = KuehnEtAl2020SInter
    MEAN_FILE = "kuehn2020/KUEHN2020_INTERFACE_GLOBAL_MEAN.csv"
    TOTAL_FILE = "kuehn2020/KUEHN2020_INTERFACE_GLOBAL_TOTAL_STDDEV.csv"
    INTER_FILE = "kuehn2020/KUEHN2020_INTERFACE_GLOBAL_INTER_EVENT_STDDEV.csv"
    INTRA_FILE = "kuehn2020/KUEHN2020_INTERFACE_GLOBAL_INTRA_EVENT_STDDEV.csv"

    def test_mean(self):
        self.check(self.MEAN_FILE, max_discrep_percentage=0.1)

    def test_std_total(self):
        self.check(self.TOTAL_FILE, max_discrep_percentage=0.1)

    def test_std_inter_event(self):
        self.check(self.INTER_FILE, max_discrep_percentage=0.1)

    def test_std_intra_event(self):
        self.check(self.INTRA_FILE, max_discrep_percentage=0.1)

    def test_mean_mb1(self):
        self.check("kuehn2020/KUEHN2020_INTERFACE_GLOBAL_Mb7.7_MEAN.csv",
                max_discrep_percentage=0.1, m_b=7.7)

    def test_mean_mb2(self):
        self.check("kuehn2020/KUEHN2020_INTERFACE_GLOBAL_Mb8.1_MEAN.csv",
                max_discrep_percentage=0.1, m_b=8.1)

    def test_mean_eps1(self):
        self.check("kuehn2020/KUEHN2020_INTERFACE_GLOBAL_epsilon1_MEAN.csv",
                max_discrep_percentage=0.1, sigma_mu_epsilon=1)

    def test_mean_eps2(self):
        self.check("kuehn2020/KUEHN2020_INTERFACE_GLOBAL_epsilon-1_MEAN.csv",
                max_discrep_percentage=0.1, sigma_mu_epsilon=-1)

# Interface Alaska
class KuehnEtAl2020SInterAlaskaTestCase(BaseGSIMTestCase):
    GSIM_CLASS = KuehnEtAl2020SInterAlaska
    MEAN_FILE = "kuehn2020/KUEHN2020_INTERFACE_ALASKA_MEAN.csv"
    TOTAL_FILE = "kuehn2020/KUEHN2020_INTERFACE_ALASKA_TOTAL_STDDEV.csv"
    INTER_FILE = "kuehn2020/KUEHN2020_INTERFACE_ALASKA_INTER_EVENT_STDDEV.csv"
    INTRA_FILE = "kuehn2020/KUEHN2020_INTERFACE_ALASKA_INTRA_EVENT_STDDEV.csv"

    def test_mean(self):
        self.check(self.MEAN_FILE, max_discrep_percentage=0.1)

    def test_std_total(self):
        self.check(self.TOTAL_FILE, max_discrep_percentage=0.1)

    def test_std_inter_event(self):
        self.check(self.INTER_FILE, max_discrep_percentage=0.1)

    def test_std_intra_event(self):
        self.check(self.INTRA_FILE, max_discrep_percentage=0.1)

# Interface Cascadia
class KuehnEtAl2020SInterCascadiaTestCase(BaseGSIMTestCase):
    GSIM_CLASS = KuehnEtAl2020SInterCascadia
    MEAN_FILE = "kuehn2020/KUEHN2020_INTERFACE_CASCADIA_MEAN.csv"
    TOTAL_FILE = "kuehn2020/KUEHN2020_INTERFACE_CASCADIA_TOTAL_STDDEV.csv"
    INTER_FILE = "kuehn2020/KUEHN2020_INTERFACE_CASCADIA_INTER_EVENT_STDDEV.csv"
    INTRA_FILE = "kuehn2020/KUEHN2020_INTERFACE_CASCADIA_INTRA_EVENT_STDDEV.csv"

    def test_mean(self):
        self.check(self.MEAN_FILE, max_discrep_percentage=0.1)

    def test_std_total(self):
        self.check(self.TOTAL_FILE, max_discrep_percentage=0.1)

    def test_std_inter_event(self):
        self.check(self.INTER_FILE, max_discrep_percentage=0.1)

    def test_std_intra_event(self):
        self.check(self.INTRA_FILE, max_discrep_percentage=0.1)

# Interface Central America and Mexico
class KuehnEtAl2020SInterCentralAmericaMexicoTestCase(BaseGSIMTestCase):
    GSIM_CLASS = KuehnEtAl2020SInterCentralAmericaMexico
    MEAN_FILE = "kuehn2020/KUEHN2020_INTERFACE_CAM_MEAN.csv"
    TOTAL_FILE = "kuehn2020/KUEHN2020_INTERFACE_CAM_TOTAL_STDDEV.csv"
    INTER_FILE = "kuehn2020/KUEHN2020_INTERFACE_CAM_INTER_EVENT_STDDEV.csv"
    INTRA_FILE = "kuehn2020/KUEHN2020_INTERFACE_CAM_INTRA_EVENT_STDDEV.csv"

    def test_mean(self):
        self.check(self.MEAN_FILE, max_discrep_percentage=0.1)

    def test_std_total(self):
        self.check(self.TOTAL_FILE, max_discrep_percentage=0.1)

    def test_std_inter_event(self):
        self.check(self.INTER_FILE, max_discrep_percentage=0.1)

    def test_std_intra_event(self):
        self.check(self.INTRA_FILE, max_discrep_percentage=0.1)

# Interface Japanclass KuehnEtAl2020SInterCascadiaTestCase(BaseGSIMTestCase):
    GSIM_CLASS = KuehnEtAl2020SInterJapan
    MEAN_FILE = "kuehn2020/KUEHN2020_INTERFACE_JAPAN_MEAN.csv"
    TOTAL_FILE = "kuehn2020/KUEHN2020_INTERFACE_JAPAN_TOTAL_STDDEV.csv"
    INTER_FILE = "kuehn2020/KUEHN2020_INTERFACE_JAPAN_INTER_EVENT_STDDEV.csv"
    INTRA_FILE = "kuehn2020/KUEHN2020_INTERFACE_JAPAN_INTRA_EVENT_STDDEV.csv"

    def test_mean(self):
        self.check(self.MEAN_FILE, max_discrep_percentage=0.1)

    def test_std_total(self):
        self.check(self.TOTAL_FILE, max_discrep_percentage=0.1)

    def test_std_inter_event(self):
        self.check(self.INTER_FILE, max_discrep_percentage=0.1)

    def test_std_intra_event(self):
        self.check(self.INTRA_FILE, max_discrep_percentage=0.1)


# Interface New Zealand
class KuehnEtAl2020SInterNewZealandTestCase(BaseGSIMTestCase):
    GSIM_CLASS = KuehnEtAl2020SInterNewZealand
    MEAN_FILE = "kuehn2020/KUEHN2020_INTERFACE_NEWZEALAND_MEAN.csv"
    TOTAL_FILE = "kuehn2020/KUEHN2020_INTERFACE_NEWZEALAND_TOTAL_STDDEV.csv"
    INTER_FILE = "kuehn2020/KUEHN2020_INTERFACE_NEWZEALAND_INTER_EVENT_STDDEV.csv"
    INTRA_FILE = "kuehn2020/KUEHN2020_INTERFACE_NEWZEALAND_INTRA_EVENT_STDDEV.csv"

    def test_mean(self):
        self.check(self.MEAN_FILE, max_discrep_percentage=0.1)

    def test_std_total(self):
        self.check(self.TOTAL_FILE, max_discrep_percentage=0.1)

    def test_std_inter_event(self):
        self.check(self.INTER_FILE, max_discrep_percentage=0.1)

    def test_std_intra_event(self):
        self.check(self.INTRA_FILE, max_discrep_percentage=0.1)

# Interface South Americaclass KuehnEtAl2020SInterCascadiaTestCase(BaseGSIMTestCase):
    GSIM_CLASS = KuehnEtAl2020SInterSouthAmerica
    MEAN_FILE = "kuehn2020/KUEHN2020_INTERFACE_SOUTHAMERICA_MEAN.csv"
    TOTAL_FILE = "kuehn2020/KUEHN2020_INTERFACE_SOUTHAMERICA_TOTAL_STDDEV.csv"
    INTER_FILE = "kuehn2020/KUEHN2020_INTERFACE_SOUTHAMERICA_INTER_EVENT_STDDEV.csv"
    INTRA_FILE = "kuehn2020/KUEHN2020_INTERFACE_SOUTHAMERICA_INTRA_EVENT_STDDEV.csv"

    def test_mean(self):
        self.check(self.MEAN_FILE, max_discrep_percentage=0.1)

    def test_std_total(self):
        self.check(self.TOTAL_FILE, max_discrep_percentage=0.1)

    def test_std_inter_event(self):
        self.check(self.INTER_FILE, max_discrep_percentage=0.1)

    def test_std_intra_event(self):
        self.check(self.INTRA_FILE, max_discrep_percentage=0.1)

# Interface Taiwan
class KuehnEtAl2020SInterTaiwanTestCase(BaseGSIMTestCase):
    GSIM_CLASS = KuehnEtAl2020SInterTaiwan
    MEAN_FILE = "kuehn2020/KUEHN2020_INTERFACE_TAIWAN_MEAN.csv"
    TOTAL_FILE = "kuehn2020/KUEHN2020_INTERFACE_TAIWAN_TOTAL_STDDEV.csv"
    INTER_FILE = "kuehn2020/KUEHN2020_INTERFACE_TAIWAN_INTER_EVENT_STDDEV.csv"
    INTRA_FILE = "kuehn2020/KUEHN2020_INTERFACE_TAIWAN_INTRA_EVENT_STDDEV.csv"

    def test_mean(self):
        self.check(self.MEAN_FILE, max_discrep_percentage=0.1)

    def test_std_total(self):
        self.check(self.TOTAL_FILE, max_discrep_percentage=0.1)

    def test_std_inter_event(self):
        self.check(self.INTER_FILE, max_discrep_percentage=0.1)

    def test_std_intra_event(self):
        self.check(self.INTRA_FILE, max_discrep_percentage=0.1)

# Interface
class KuehnEtAl2020SSlabTestCase(BaseGSIMTestCase):
    GSIM_CLASS = KuehnEtAl2020SSlab
    MEAN_FILE = "kuehn2020/KUEHN2020_INSLAB_GLOBAL_MEAN.csv"
    TOTAL_FILE = "kuehn2020/KUEHN2020_INSLAB_GLOBAL_TOTAL_STDDEV.csv"
    INTER_FILE = "kuehn2020/KUEHN2020_INSLAB_GLOBAL_INTER_EVENT_STDDEV.csv"
    INTRA_FILE = "kuehn2020/KUEHN2020_INSLAB_GLOBAL_INTRA_EVENT_STDDEV.csv"

    def test_mean(self):
        self.check(self.MEAN_FILE, max_discrep_percentage=0.1)

    def test_std_total(self):
        self.check(self.TOTAL_FILE, max_discrep_percentage=0.1)

    def test_std_inter_event(self):
        self.check(self.INTER_FILE, max_discrep_percentage=0.1)

    def test_std_intra_event(self):
        self.check(self.INTRA_FILE, max_discrep_percentage=0.1)

    def test_mean_mb1(self):
        self.check("kuehn2020/KUEHN2020_INSLAB_GLOBAL_Mb7.4_MEAN.csv",
                max_discrep_percentage=0.1, m_b=7.4)

    def test_mean_mb2(self):
        self.check("kuehn2020/KUEHN2020_INSLAB_GLOBAL_Mb7.8_MEAN.csv",
                max_discrep_percentage=0.1, m_b=7.8)

    def test_mean_eps1(self):
        self.check("kuehn2020/KUEHN2020_INSLAB_GLOBAL_epsilon1_MEAN.csv",
                max_discrep_percentage=0.1, sigma_mu_epsilon=1)

    def test_mean_eps2(self):
        self.check("kuehn2020/KUEHN2020_INSLAB_GLOBAL_epsilon-1_MEAN.csv",
                max_discrep_percentage=0.1, sigma_mu_epsilon=-1)

# Intraslab Alaska
class KuehnEtAl2020SSlabAlaskaTestCase(BaseGSIMTestCase):
    GSIM_CLASS = KuehnEtAl2020SSlabAlaska
    MEAN_FILE = "kuehn2020/KUEHN2020_INSLAB_ALASKA_MEAN.csv"
    TOTAL_FILE = "kuehn2020/KUEHN2020_INSLAB_ALASKA_TOTAL_STDDEV.csv"
    INTER_FILE = "kuehn2020/KUEHN2020_INSLAB_ALASKA_INTER_EVENT_STDDEV.csv"
    INTRA_FILE = "kuehn2020/KUEHN2020_INSLAB_ALASKA_INTRA_EVENT_STDDEV.csv"

    def test_mean(self):
        self.check(self.MEAN_FILE, max_discrep_percentage=0.1)

    def test_std_total(self):
        self.check(self.TOTAL_FILE, max_discrep_percentage=0.1)

    def test_std_inter_event(self):
        self.check(self.INTER_FILE, max_discrep_percentage=0.1)

    def test_std_intra_event(self):
        self.check(self.INTRA_FILE, max_discrep_percentage=0.1)

# Intraslab Cascadia
class KuehnEtAl2020SSlabCascadiaTestCase(BaseGSIMTestCase):
    GSIM_CLASS = KuehnEtAl2020SSlabCascadia
    MEAN_FILE = "kuehn2020/KUEHN2020_INSLAB_CASCADIA_MEAN.csv"
    TOTAL_FILE = "kuehn2020/KUEHN2020_INSLAB_CASCADIA_TOTAL_STDDEV.csv"
    INTER_FILE = "kuehn2020/KUEHN2020_INSLAB_CASCADIA_INTER_EVENT_STDDEV.csv"
    INTRA_FILE = "kuehn2020/KUEHN2020_INSLAB_CASCADIA_INTRA_EVENT_STDDEV.csv"

    def test_mean(self):
        self.check(self.MEAN_FILE, max_discrep_percentage=0.1)

    def test_std_total(self):
        self.check(self.TOTAL_FILE, max_discrep_percentage=0.1)

    def test_std_inter_event(self):
        self.check(self.INTER_FILE, max_discrep_percentage=0.1)

    def test_std_intra_event(self):
        self.check(self.INTRA_FILE, max_discrep_percentage=0.1)

# Interface Central America and Mexico
class KuehnEtAl2020SInterCentralAmericaMexicoTestCase(BaseGSIMTestCase):
    GSIM_CLASS = KuehnEtAl2020SSlabCentralAmericaMexico
    MEAN_FILE = "kuehn2020/KUEHN2020_INSLAB_CAM_MEAN.csv"
    TOTAL_FILE = "kuehn2020/KUEHN2020_INSLAB_CAM_TOTAL_STDDEV.csv"
    INTER_FILE = "kuehn2020/KUEHN2020_INSLAB_CAM_INTER_EVENT_STDDEV.csv"
    INTRA_FILE = "kuehn2020/KUEHN2020_INSLAB_CAM_INTRA_EVENT_STDDEV.csv"

    def test_mean(self):
        self.check(self.MEAN_FILE, max_discrep_percentage=0.1)

    def test_std_total(self):
        self.check(self.TOTAL_FILE, max_discrep_percentage=0.1)

    def test_std_inter_event(self):
        self.check(self.INTER_FILE, max_discrep_percentage=0.1)

    def test_std_intra_event(self):
        self.check(self.INTRA_FILE, max_discrep_percentage=0.1)

# Interface Japanclass KuehnEtAl2020SInterCascadiaTestCase(BaseGSIMTestCase):
    GSIM_CLASS = KuehnEtAl2020SSlabJapan
    MEAN_FILE = "kuehn2020/KUEHN2020_INSLAB_JAPAN_MEAN.csv"
    TOTAL_FILE = "kuehn2020/KUEHN2020_INSLAB_JAPAN_TOTAL_STDDEV.csv"
    INTER_FILE = "kuehn2020/KUEHN2020_INSLAB_JAPAN_INTER_EVENT_STDDEV.csv"
    INTRA_FILE = "kuehn2020/KUEHN2020_INSLAB_JAPAN_INTRA_EVENT_STDDEV.csv"

    def test_mean(self):
        self.check(self.MEAN_FILE, max_discrep_percentage=0.1)

    def test_std_total(self):
        self.check(self.TOTAL_FILE, max_discrep_percentage=0.1)

    def test_std_inter_event(self):
        self.check(self.INTER_FILE, max_discrep_percentage=0.1)

    def test_std_intra_event(self):
        self.check(self.INTRA_FILE, max_discrep_percentage=0.1)


# Intraslab New Zealand
class KuehnEtAl2020SSlabNewZealandTestCase(BaseGSIMTestCase):
    GSIM_CLASS = KuehnEtAl2020SSlabNewZealand
    MEAN_FILE = "kuehn2020/KUEHN2020_INSLAB_NEWZEALAND_MEAN.csv"
    TOTAL_FILE = "kuehn2020/KUEHN2020_INSLAB_NEWZEALAND_TOTAL_STDDEV.csv"
    INTER_FILE = "kuehn2020/KUEHN2020_INSLAB_NEWZEALAND_INTER_EVENT_STDDEV.csv"
    INTRA_FILE = "kuehn2020/KUEHN2020_INSLAB_NEWZEALAND_INTRA_EVENT_STDDEV.csv"

    def test_mean(self):
        self.check(self.MEAN_FILE, max_discrep_percentage=0.1)

    def test_std_total(self):
        self.check(self.TOTAL_FILE, max_discrep_percentage=0.1)

    def test_std_inter_event(self):
        self.check(self.INTER_FILE, max_discrep_percentage=0.1)

    def test_std_intra_event(self):
        self.check(self.INTRA_FILE, max_discrep_percentage=0.1)

# Intraslab South America
class KuehnEtAl2020SSlabSouthAmericaTestCase(BaseGSIMTestCase):
    GSIM_CLASS = KuehnEtAl2020SSlabSouthAmerica
    MEAN_FILE = "kuehn2020/KUEHN2020_INSLAB_SOUTHAMERICA_MEAN.csv"
    TOTAL_FILE = "kuehn2020/KUEHN2020_INSLAB_SOUTHAMERICA_TOTAL_STDDEV.csv"
    INTER_FILE = "kuehn2020/KUEHN2020_INSLAB_SOUTHAMERICA_INTER_EVENT_STDDEV.csv"
    INTRA_FILE = "kuehn2020/KUEHN2020_INSLAB_SOUTHAMERICA_INTRA_EVENT_STDDEV.csv"

    def test_mean(self):
        self.check(self.MEAN_FILE, max_discrep_percentage=0.1)

    def test_std_total(self):
        self.check(self.TOTAL_FILE, max_discrep_percentage=0.1)

    def test_std_inter_event(self):
        self.check(self.INTER_FILE, max_discrep_percentage=0.1)

    def test_std_intra_event(self):
        self.check(self.INTRA_FILE, max_discrep_percentage=0.1)

# Intraslab Taiwan
class KuehnEtAl2020SSlabTaiwanTestCase(BaseGSIMTestCase):
    GSIM_CLASS = KuehnEtAl2020SSlabTaiwan
    MEAN_FILE = "kuehn2020/KUEHN2020_INSLAB_TAIWAN_MEAN.csv"
    TOTAL_FILE = "kuehn2020/KUEHN2020_INSLAB_TAIWAN_TOTAL_STDDEV.csv"
    INTER_FILE = "kuehn2020/KUEHN2020_INSLAB_TAIWAN_INTER_EVENT_STDDEV.csv"
    INTRA_FILE = "kuehn2020/KUEHN2020_INSLAB_TAIWAN_INTRA_EVENT_STDDEV.csv"

    def test_mean(self):
        self.check(self.MEAN_FILE, max_discrep_percentage=0.1)

    def test_std_total(self):
        self.check(self.TOTAL_FILE, max_discrep_percentage=0.1)

    def test_std_inter_event(self):
        self.check(self.INTER_FILE, max_discrep_percentage=0.1)

    def test_std_intra_event(self):
        self.check(self.INTRA_FILE, max_discrep_percentage=0.1)

