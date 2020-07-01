# Code and data associated with "Reductions in commuting mobility predict geographic differences in SARS-CoV-2 prevalence in New York City"

## Stephen M. Kissler, Nishant Kishore, Malavika Prabhu, Dena Goffman, Yaakov Beilin, Ruth Landau, Cynthia Gyamfi-Bannerman, Brian T. Bateman, Jon Snyder, Armin S. Razavi, Daniel Katz, Jonathan Gal, Angela Bianco, Joanne Stone, Daniel Larremore, Caroline O. Buckee, Yonatan H. Grad

https://dash.harvard.edu/bitstream/handle/1/42665370/Kissler_etal_NYC_mobility.pdf?sequence=1&isAllowed=y

Analysis was run using R version 3.6.2 on a 2017 MacBook Pro using an 2.3 GHz Intel Core i5 processor.

__figuremaker.R:__ Main file for reproducing the figures and analysis from the article. Expected runtime: 5 min.

__serofunctions2.R:__ Code for inferring population prevalence given test sensitivity and specificity; see https://github.com/LarremoreLab/covid_serological_sampling and https://larremorelab.github.io/covid19testgroup

__borough_zip.csv:__ Mapping from boroughs to ZIP codes

__borough_zip3.csv:__ Mapping from boroughs to 3-digit ZIP code prefixes

__borough_zip3.csv:__ Mapping from boroughs to 3-digit ZIP code prefixes, separating North and South Queens

__nyc_hospital_data.csv:__ Data from the New York City hospitals

All code is available under the GNU General Public License, version 3, included in this repository under the following terms: 
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
