import re
import sys
import requests
import time
import random
import pandas as pd
from bs4 import BeautifulSoup
from tqdm import tqdm

def find_re(root, expr):
  node = root.find(text=expr)
  return expr.match(node.text) if node else None

def shuffle(L):
  random.shuffle(L)
  return L

base = 'https://www.drugs.com'
backoff = 0.05
au_pregnancy_cat = re.compile(r'^AU TGA pregnancy category: (.+)$', re.IGNORECASE)
us_pregnancy_cat = re.compile(r'^US FDA pregnancy category: (.+)$', re.IGNORECASE)
us_pregnancy_cat_2 = re.compile(r'.+pregnancy category (\w) by the FDA.+', re.IGNORECASE)

# Pregnancy home page
html = requests.get('https://www.drugs.com/pregnancy/').text
#  Get all alphabetical subsections page links
a_z_sections = [f"{base}{a.get('href')}" for a in BeautifulSoup(html).find(text='Browse A-Z').parent.parent.find_all('a')]
drug_pages = {}

for a_z_section in tqdm(a_z_sections):
  time.sleep(backoff)
  html = requests.get(a_z_section).text
  root = BeautifulSoup(html)
  for a in root.find(text='Alphabetical').parent.find_next_sibling('ul').find_all('a'):
    if a.get('href').startswith('/pregnancy/'):
      drug_pages[f"{base}{a.get('href')}"] = a.text

drug_pregnancy_cat = {}

for drug_page, drug in tqdm([
  (drug_page, drug)
  for drug_page, drug in drug_pages.items()
  if drug_page not in drug_pregnancy_cat
]):
  time.sleep(backoff)
  html = requests.get(drug_page).text
  root = BeautifulSoup(html)
  au_m = find_re(root, au_pregnancy_cat)
  au = au_m.group(1) if au_m else None
  us_m = find_re(root, us_pregnancy_cat)
  us = us_m.group(1) if us_m else None
  if us is None:
    us_m = find_re(root, us_pregnancy_cat_2)
    us = us_m.group(1) if us_m else None
  drug_pregnancy_cat[drug_page] = dict(drug=drug, au=au, us=us)

df = pd.DataFrame.from_dict(drug_pregnancy_cat, orient='index')

df['au-cleaned'] = df.au.str.strip(r'[ \.]').str.upper().replace({
  'EXEMPT': None,
  'NOT FORMALLY ASSIGNED TO A PREGNANCY CATEGORY': None,
  'NOT FORMALLY ASSIGNED TO A PREGNANCY CATEGORY.  THIS CLASS OF DRUGS IS GENERALLY EXEMPT FROM PREGNANCY CLASSIFICATION': None,
  'NOT ASSIGNED': None,
  'NOT ASSIGNED; THIS CLASS OF DRUGS IS GENERALLY EXEMPT FROM PREGNANCY CLASSIFICATION': None,
  'B3 (ESTRADIOL); CATEGORY A (PROGESTERONE)': 'B3',
  'A (SHORT-TERM THERAPY)': 'A',
  'CATEGORY D': 'D',
  'A (ORAL, RECTAL FOAM); C (PARENTERAL)': 'A',
  'THIS DRUG IS NOT INDICATED FOR USE IN FEMALES AND HAS NOT BEEN ASSIGNED TO A TGA PREGNANCY CATEGORY': None,
  'GEL/OINTMENT FORMULATIONS: B1; FOAM FORMULATION: B3': 'B3',
  'B 3': 'B3',
  'CATEGORY B1': 'B1',
  'A (ORAL); C (PARENTERAL)': 'C',
  'A (INHALATION); B3 (ORAL)': 'B3',
  'CATEGORY A (ACETAMINOPHEN/PARACETAMOL)': 'A',
  'NOT AVAILABLE': None,
})
df['au-cleaned'].value_counts()
df['us-cleaned'] = df.us.str.strip(r'[ \.]').str.upper().replace({
  'CATEGORY X': 'X',
  'NOT ASSIGNED': None,
  'NOT AVAILABLE': None,
  'NOT FORMALLY ASSIGNED TO A PREGNANCY CATEGORY': None,
  'NOT FORMALLY ASSIGNED A PREGNANCY CATEGORY': None,
  'NOT FORMALLY ASSIGNED TO A PREGNANCY CATEGORY.  NOT FORMALLY ASSIGNED TO A PREGNANCY CATEGORY': None,
  'NOT FORMALLY ASSIGNED TO A PREGNANCY CATEGORY DUE TO PRODUCT BEING AVAILABLE OVER-THE-COUNTER.  OVER-THE-COUNTER WARNING: IT IS ESPECIALLY IMPORTANT NOT TO USE ASPIRIN DURING THE LAST 3 MONTHS OF PREGNANCY UNLESS DEFINITELY DIRECTED TO DO SO BY A DOCTOR BECAUSE IT MAY CAUSE PROBLEMS IN THE UNBORN CHILD OR COMPLICATIONS DURING DELIVERY': None,
  'NOT FORMALLY ASSIGNED': None,
  'UNASSIGNED': None,
  'C PRIOR TO 30 WEEKS GESTATION': 'C',
  'C (IMMEDIATE-RELEASE); NOT ASSIGNED (EXTENDED-RELEASE)': 'C',
  'CONTRAINDICATION': None,
  'D (IV); NOT ASSIGNED (TABLETS)': 'D',
  'C/D': 'D',
  'X (PSORIASIS AND RHEUMATOID ARTHRITIS); NOT ASSIGNED (ALL OTHER CONDITIONS)': 'X',
  'X (FOR DRUG-INDUCED METHEMOGLOBINEMIA); NOT ASSIGNED (FOR ACQUIRED METHEMOGLOBINEMIA)': 'X',
  'X (BRAND KORLYM)': 'X',
  'C (NITROPRESS); NOT ASSIGNED': 'C',
  'D/B (MANUFACTURER DEPENDENT)': 'D',
  'NOT ASSIGNED (EXTENDED-RELEASE)': None,
  'C (1% CREAM)': 'C',
  'C (TABLETS); NOT ASSIGNED (INJECTABLE SOLUTION)': 'C',
  'D (THIRD TRIMESTER); C (FIRST AND SECOND TRIMESTER)': 'D',
  'C (OPHTHALMIC OINTMENT)': 'C',
  'N': None,
  'C (0.45% AND 0.5% OPHTHALMIC SOLUTIONS); NOT ASSIGNED (0.4% OPHTHALMIC SOLUTION)': 'C',
  'D IN PATIENTS WITH ADVANCED BREAST CANCER; X IN PATIENTS WITH ENDOMETRIOSIS AND ENDOMETRIAL THINNING': 'X',
  'C/B (MANUFACTURER DEPENDENT)': 'C',
  'NOT ASSIGNED (U-100); CATEGORY B (U-500)': 'B',
  'C (PRIOR TO 30 WEEKS GESTATION)': 'C',
  'D (ALL INDICATIONS EXCEPT VAGINAL CANDIDIASIS)': 'D',
  'B (IV); NOT ASSIGNED (ORAL)': 'B',
  'B (3% GEL); C PRIOR TO 30 WEEKS GESTATION; D STARTING AT 30 WEEKS GESTATION (TOPICAL SOLUTIONS, 1% GEL)': 'D',
  'C (TABLETS)': 'C',
  'C (ORAL); NOT ASSIGNED (INJECTABLE)': 'C',
  'C (PERIODONTAL CHIP)': 'C',
  'D (SECOND AND THIRD TRIMESTERS); C (FIRST TRIMESTER)': 'D',
  'IV: NOT ASSIGNED; ORAL: D': 'D',
  'C (0.07% AND 0.09% OPHTHALMIC SOLUTIONS); NOT ASSIGNED (0.075% OPHTHALMIC SOLUTION)': 'C',
  'C (SUPPOSITORIES)': 'C',
  'C (TABLETS); B (ORAL INHALATION); NOT ASSIGNED (PARENTERAL)': 'C',
})
df['us-cleaned'].value_counts()

if __name__ == '__main__':
  assert sys.argv[1].endswith('.tsv')
  df.rename_axis('url').reset_index().set_index('drug')[[
    'us-cleaned', 'au-cleaned', 'us', 'au', 'url'
  ]].to_csv(sys.argv[1], sep='\t')
