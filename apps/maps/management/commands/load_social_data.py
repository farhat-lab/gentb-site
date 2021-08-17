"""
Load the social determinant of health data from the World Health Organization and the World Bank
"""
import csv
import requests
import sys
import os
from urllib.request import urlopen
import io
from zipfile import ZipFile
from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
sys.path.insert(0, '.')
sys.path.insert(0, '..')

try:
    import manage # pylint: disable=unused-import
except ImportError as err:
    sys.stderr.write("Could not run script! Is manage.py not in the current"\
        "working directory, or is the environment not configured?:\n"\
        "{}\n".format(err))
    sys.exit(1)

from apps.mutations.utils import csv_merge
from apps.maps.models import Country, CountryHealth
from apps.mutations.utils import StatusBar

class Command(BaseCommand):
    def handle(self, verbose=True, **options):
        burden_url = 'https://extranet.who.int/tme/generateCSV.asp?ds=estimates'
        mdr_url = "https://extranet.who.int/tme/generateCSV.asp?ds=mdr_rr_estimates"
        budget_url = "https://extranet.who.int/tme/generateCSV.asp?ds=budget"
        latent_url = "https://extranet.who.int/tme/generateCSV.asp?ds=ltbi_estimates"
        pop_dens_url = "http://api.worldbank.org/v2/en/indicator/EN.POP.DNST?downloadformat=csv"
        wealth_url = "http://api.worldbank.org/v2/en/indicator/NW.TOW.PC?downloadformat=csv"
        gdp_url = "http://api.worldbank.org/v2/en/indicator/NY.GDP.MKTP.CD?downloadformat=csv"
        avika_data = 'Supp Table 3 who_corrected_snp10_6.4.21.csv'


        with requests.Session() as s:
            #load tb burden data for hiv coincidence and overall tb incidence
            download_burden = s.get(burden_url)
            decoded_content_burden = download_burden.content.decode('utf-8')
            cr_burden = csv.reader(decoded_content_burden.splitlines(), delimiter=',')
            csv_list_burden = list(cr_burden)
            header_burden = csv_list_burden[0]
            country_burden_index = header_burden.index("iso3")
            year_burden_index = header_burden.index("year")
            hiv_index = header_burden.index('e_inc_tbhiv_100k')
            incidence_index = header_burden.index("e_inc_100k")
            for row in csv_list_burden[1:]:
                if row[year_burden_index]=='2018':
                    try:
                        country = Country.objects.get(iso3=row[country_burden_index])
                        health = country.health
                    except CountryHealth.DoesNotExist:
                        health = CountryHealth(country=country)
                    except Country.DoesNotExist:
                        print("Can't find country {}".format(row[country_burden_index]))
                        continue
                    health.hiv_incidence2018 = row[hiv_index]
                    health.all_tb_incidence2018 = row[incidence_index]
                    health.save()
                    print("Saved WHO TB burden data for {}".format(country))

        with requests.Session() as s:
            #load who data for estimated mdr
            download_mdr = s.get(mdr_url)
            decoded_content_mdr = download_mdr.content.decode('utf-8')
            cr_mdr = csv.reader(decoded_content_mdr.splitlines(), delimiter=',')
            csv_list_mdr = list(cr_mdr)
            header_mdr = csv_list_mdr[0]
            country_mdr_index = header_mdr.index("iso3")
            est_mdr_index = header_mdr.index('e_rr_pct_new')
            for row in csv_list_mdr[1:]:
                try:
                    country = Country.objects.get(iso3=row[country_mdr_index])
                    health = country.health
                except CountryHealth.DoesNotExist:
                    health = CountryHealth(country=country)
                except Country.DoesNotExist:
                    print("Can't find country {}".format(row[country_mdr_index]))
                    continue
                health.est_mdr = row[est_mdr_index]
                health.save()
                print("Saved WHO est mdr data for {}".format(country))

        with requests.Session() as s:
            #load who budget data for total funding
            download_budget = s.get(budget_url)
            decoded_content_budget = download_budget.content.decode('utf-8')
            cr_budget = csv.reader(decoded_content_budget.splitlines(), delimiter=',')
            csv_list_budget = list(cr_budget)
            header_budget = csv_list_budget[0]
            country_budget_index = header_budget.index("iso3")
            country_funding_index = header_budget.index("cf_tot_sources")
            for row in csv_list_budget[1:]:
                try:
                    country = Country.objects.get(iso3=row[country_budget_index])
                    health = country.health
                except CountryHealth.DoesNotExist:
                    health = CountryHealth(country=country)
                except Country.DoesNotExist:
                    print("Can't find country {}".format(row[country_budget_index]))
                    continue
                health.total_funding = row[country_funding_index]
                health.save()
                print("Saved WHO TB budget data for {}".format(country))

        with requests.Session() as s:
            # load who tb latent data for estimated household size
            download_latent = s.get(latent_url)
            decoded_content_latent = download_latent.content.decode('utf-8')
            cr_latent = csv.reader(decoded_content_latent.splitlines(), delimiter=',')
            csv_list_latent = list(cr_latent)
            header_latent = csv_list_latent[0]
            country_latent_index = header_latent.index("iso3")
            household_index = header_latent.index("e_hh_size")
            for row in csv_list_latent[1:]:
                try:
                    country = Country.objects.get(iso3=row[country_latent_index])
                    health = country.health
                except CountryHealth.DoesNotExist:
                    health = CountryHealth(country=country)
                except Country.DoesNotExist:
                    print("Can't find country {}".format(row[country_latent_index]))
                    continue
                health.household = row[household_index]
                health.save()
                print("Saved WHO TB latent data for {}".format(country))

        # load world bank population density data
        pop_package = io.BytesIO(urlopen(pop_dens_url).read())
        pop_z = ZipFile(pop_package, 'r')
        DATA_DIR = os.path.join(settings.DATA_ROOT, 'maps')
        pop_pwd = os.path.join(DATA_DIR, pop_z.namelist()[1])
        with open(pop_pwd, 'w') as fhl:
            fhl.write(pop_z.read(pop_z.namelist()[1]).decode(('utf-8')))
        with open(pop_pwd, 'r') as pop_csv:
            rows = csv.reader(pop_csv)
            next(rows)
            next(rows)
            next(rows)
            next(rows)
            header = next(rows)
            country_index = header.index("Country Code")
            year_index = header.index("2018")
            for row in rows:
                try:
                    country = Country.objects.get(iso3=row[country_index])
                    health = country.health
                except CountryHealth.DoesNotExist:
                    health = CountryHealth(country=country)
                except Country.DoesNotExist:
                    print("Can't find country {}".format(row[country_index]))
                    continue
                health.pop_dens = row[year_index]
                health.save()
                print("Saved World Bank population density data for {}".format(country))

        # load world bank total wealth per capita data
        wealth_package = io.BytesIO(urlopen(wealth_url).read())
        wealth_z = ZipFile(wealth_package, 'r')
        DATA_DIR = os.path.join(settings.DATA_ROOT, 'maps')
        wealth_pwd = os.path.join(DATA_DIR, wealth_z.namelist()[1])
        with open(wealth_pwd, 'w') as fhl:
            fhl.write(wealth_z.read(wealth_z.namelist()[1]).decode(('utf-8')))
        with open(wealth_pwd, 'r') as wealth_csv:
            rows = csv.reader(wealth_csv)
            next(rows)
            next(rows)
            next(rows)
            next(rows)
            header = next(rows)
            country_index = header.index("Country Code")
            year_index = header.index("2014")
            for row in rows:
                try:
                    country = Country.objects.get(iso3=row[country_index])
                    health = country.health
                except CountryHealth.DoesNotExist:
                    health = CountryHealth(country=country)
                except Country.DoesNotExist:
                    print("Can't find country {}".format(row[country_index]))
                    continue
                health.total_wealth = row[year_index]
                health.save()
                print("Saved World Bank total wealth data for {}".format(country))

        # load world bank gdp data
        gdp_package = io.BytesIO(urlopen(gdp_url).read())
        gdp_z = ZipFile(gdp_package, 'r')
        DATA_DIR = os.path.join(settings.DATA_ROOT, 'maps')
        gdp_pwd = os.path.join(DATA_DIR, gdp_z.namelist()[1])
        with open(gdp_pwd, 'w') as fhl:
            fhl.write(gdp_z.read(gdp_z.namelist()[1]).decode(('utf-8')))
        with open(gdp_pwd, 'r') as gdp_csv:
            rows = csv.reader(gdp_csv)
            next(rows)
            next(rows)
            next(rows)
            next(rows)
            header = next(rows)
            country_index = header.index("Country Code")
            year_index = header.index("2019")
            for row in rows:
                try:
                    country = Country.objects.get(iso3=row[country_index])
                    health = country.health
                except CountryHealth.DoesNotExist:
                    health = CountryHealth(country=country)
                except Country.DoesNotExist:
                    print("Can't find country {}".format(row[country_index]))
                health.world_bank_gdp = row[year_index]
                health.save()
                print("Saved World Bank GDP data for {}".format(country))

        # load Avika's data
        DATA_DIR = os.path.join(settings.DATA_ROOT, 'maps')
        avika_pwd = os.path.join(DATA_DIR, avika_data)
        with open(avika_pwd, 'r') as avika_csv:
            rows = csv.reader(avika_csv)
            header = next(rows)
            country_index = header.index("country")
            drug_index = header.index("drug")
            """rr_drug = {
                    "amk" : health.amk_rr,
                    "cap" : health.cap_rr,
                    "cip" : health.cip_rr,
                    "emb" : health.emb_rr,
                    "eth" : health.eth_rr,
                    "inh" : health.inh_rr,
                    "kan" : health.kan_rr,
                    "levo" : health.levo_rr,
                    "oflx" : health.oflx_rr,
                    "pas" : health.pas_rr,
                    "pza" : health.pza_rr,
                    "rif" : health.rif_rr,
                    "str" : health.str_rr 
             }"""
            for row in rows:
                try:
                    country = Country.objects.get(name=row[country_index])
                    health = country.health
                except CountryHealth.DoesNotExist:
                    health = CountryHealth(country=country)
                except Country.DoesNotExist:
                    print("Can't find country {}".format(row[country_index]))
                    continue

                if (row[drug_index] == "amk") : 
                    health.amk_resistant = row[header.index("gentb_snp10_res")]
                    health.amk_rr = row[header.index("gentb_snp10_res_rr")]
                    health.amk_rs = row[header.index("gentb_snp10_res_rs")]
                    health.amk_mean_rr = row[header.index("mean_snp10_rr")]
                #rr = row[header.index("p_a_r")]
                    health.amk_mean_rs = row[header.index("mean_snp10_rs")]

                if (row[drug_index] == "cap") :
                    health.cap_resistant = row[header.index("gentb_snp10_res")]
                    health.cap_rr = row[header.index("gentb_snp10_res_rr")]
                    health.cap_rs = row[header.index("gentb_snp10_res_rs")]
                    health.cap_mean_rr = row[header.index("mean_snp10_rr")]
                #rr = row[header.index("p_a_r")]
                    health.cap_mean_rs = row[header.index("mean_snp10_rs")]

                if (row[drug_index] == "cip") :
                    health.cip_resistant = row[header.index("gentb_snp10_res")]
                    health.cip_rr = row[header.index("gentb_snp10_res_rr")]
                    health.cip_rs = row[header.index("gentb_snp10_res_rs")]
                    health.cip_mean_rr = row[header.index("mean_snp10_rr")]
                #rr = row[header.index("p_a_r")]
                    health.cip_mean_rs = row[header.index("mean_snp10_rs")]

                if (row[drug_index] == "emb") :
                    health.emb_resistant = row[header.index("gentb_snp10_res")]
                    health.emb_rr = row[header.index("gentb_snp10_res_rr")]
                    health.emb_rs = row[header.index("gentb_snp10_res_rs")]
                    health.emb_mean_rr = row[header.index("mean_snp10_rr")]
                #rr = row[header.index("p_a_r")]
                    health.emb_mean_rs = row[header.index("mean_snp10_rs")]

                if (row[drug_index] == "eth") :
                    health.eth_resistant = row[header.index("gentb_snp10_res")]
                    health.eth_rr = row[header.index("gentb_snp10_res_rr")]
                    health.eth_rs = row[header.index("gentb_snp10_res_rs")]
                    health.eth_mean_rr = row[header.index("mean_snp10_rr")]
                #rr = row[header.index("p_a_r")]
                    health.eth_mean_rs = row[header.index("mean_snp10_rs")]

                if (row[drug_index] == "inh") :
                    health.inh_resistant = row[header.index("gentb_snp10_res")]
                    health.inh_rr = row[header.index("gentb_snp10_res_rr")]
                    health.inh_rs = row[header.index("gentb_snp10_res_rs")]
                    health.inh_mean_rr = row[header.index("mean_snp10_rr")]
                #rr = row[header.index("p_a_r")]
                    health.inh_mean_rs = row[header.index("mean_snp10_rs")]

                if (row[drug_index] == "kan") :
                    health.kan_resistant = row[header.index("gentb_snp10_res")]
                    health.kan_rr = row[header.index("gentb_snp10_res_rr")]
                    health.kan_rs = row[header.index("gentb_snp10_res_rs")]
                    health.kan_mean_rr = row[header.index("mean_snp10_rr")]
                #rr = row[header.index("p_a_r")]
                    health.kan_mean_rs = row[header.index("mean_snp10_rs")]

                if (row[drug_index] == "levo") :
                    health.levo_resistant = row[header.index("gentb_snp10_res")]
                    health.levo_rr = row[header.index("gentb_snp10_res_rr")]
                    health.levo_rs = row[header.index("gentb_snp10_res_rs")]
                    health.levo_mean_rr = row[header.index("mean_snp10_rr")]
                #rr = row[header.index("p_a_r")]
                    health.levo_mean_rs = row[header.index("mean_snp10_rs")]

                if (row[drug_index] == "oflx") :
                    health.oflx_resistant = row[header.index("gentb_snp10_res")]
                    health.oflx_rr = row[header.index("gentb_snp10_res_rr")]
                    health.oflx_rs = row[header.index("gentb_snp10_res_rs")]
                    health.oflx_mean_rr = row[header.index("mean_snp10_rr")]
                #rr = row[header.index("p_a_r")]
                    health.oflx_mean_rs = row[header.index("mean_snp10_rs")]

                if (row[drug_index] == "pas") :
                    health.pas_resistant = row[header.index("gentb_snp10_res")]
                    health.pas_rr = row[header.index("gentb_snp10_res_rr")]
                    health.pas_rs = row[header.index("gentb_snp10_res_rs")]
                    health.pas_mean_rr = row[header.index("mean_snp10_rr")]
                #rr = row[header.index("p_a_r")]
                    health.pas_mean_rs = row[header.index("mean_snp10_rs")]

                if (row[drug_index] == "pza") :
                    health.pza_resistant = row[header.index("gentb_snp10_res")]
                    health.pza_rr = row[header.index("gentb_snp10_res_rr")]
                    health.pza_rs = row[header.index("gentb_snp10_res_rs")]
                    health.pza_mean_rr = row[header.index("mean_snp10_rr")]
                #rr = row[header.index("p_a_r")]
                    health.pza_mean_rs = row[header.index("mean_snp10_rs")]

                if (row[drug_index] == "rif") :
                    health.rif_resistant = row[header.index("gentb_snp10_res")]
                    health.rif_rr = row[header.index("gentb_snp10_res_rr")]
                    health.rif_rs = row[header.index("gentb_snp10_res_rs")]
                    health.rif_mean_rr = row[header.index("mean_snp10_rr")]
                #rr = row[header.index("p_a_r")]
                    health.rif_mean_rs = row[header.index("mean_snp10_rs")]

                if (row[drug_index] == "str") : 
                    health.str_resistant = row[header.index("gentb_snp10_res")]
                    health.str_rr = row[header.index("gentb_snp10_res_rr")]
                    health.str_rs = row[header.index("gentb_snp10_res_rs")]
                    health.str_mean_rr = row[header.index("mean_snp10_rr")]
                #rr = row[header.index("p_a_r")]
                    health.str_mean_rs = row[header.index("mean_snp10_rs")]
                """
                drugObj = lookup.get(row[drug_index],None)
                print(row[drug_index])
                print(drugObj)
                drugObj.name = row[drug_index]
                drugObj.resistant = row[header.index("gentb_snp10_res")]
                drugObj.rr = row[header.index("gentb_snp10_res_rr")]
                drugObj.rs = row[header.index("gentb_snp10_res_rs")]
                drugObj.mean_rr = row[header.index("mean_snp10_rr")]
                drugObj.rr = row[header.index("p_a_r")]
                drugObj.mean_rs = row[header.index("mean_snp10_rs")]
                """
                #health.percent_resistance = row[header.index('p_a_r')]
                #health.confidence_interval = row[header.index('var_a_r')]
                health.sample_size = row[header.index('gentb_snp10_n')]
                health.resistant_isolates = row[header.index("gentb_snp10_n_rr")]
                health.suseptable_isolates = row[header.index("gentb_snp10_n_rs")]
                #print(row[header.index('gentb_snp10_n')])
             # health.isonized_resistance = row[header.index('mean_snp10_rs')]
                health.save()
                print("Saved Avika Data for {}".format(country))
                print("drug {}".format(row[drug_index]))

