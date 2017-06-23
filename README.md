# REU-Gurobi
Gurobi/Python documents for 2017 REU program

Reading From Excel/Data Formatting
The best libraries for this are pandas or openpyxl. I mostly use pandas, but if the ability to specify a specific cell/cell range is important, openpyxl may be easier. 

Reads and parses in the file into a pandas DataFrame object:
file = pandas.ExcelFile('.\PATH\data.xlsx')
my_dataframe = file.parse(<args>)

Or, to read and parse in a single command:
my_dataframe = pandas.read_excel('.\PATH\data.xlsx', <args>)

See documentation for pandas.read_excel for an explanation of parse arguments. 

Gurobi doesn't recognize DataFrames, so they need to be converted to dictionaries or lists.
Also, strings are read in unicode. I think recognizes unicode, so this is partly up to personal preference, but you need to be consistent so that dictionary keys match.

For a single column you can just use DataFrame.values.flatten()
Example 1:
  hs_df = file.parse(header=1,skipfooter=677,parse_cols='M')
  u_hs = hs_df.values.flatten()
  hs = [str(x) for x in u_hs] # <-- Converts keys from unicode to string.

For more than one column, use DataFrame.to_dict(). See the documentation for more flexibility.
The details will change depending on how the Excel sheet is structured, how many variables you have, and what data structure you want to end up with, but other DataFrames also have other manipulation methods that may be helpful, such as set_index().

For an Excel Sheet with the structure
Sites   Distance
HS1     100
HS2     200
..      ...
to_dict can be called directly once the index is correct.
Exmaple 2:
  u_tau = tau_df.set_index('Sites').to_dict()
  tau = {str(k): v for k, v in u_tau.items()} # <-- dict comprehension is used to convert the keys from unicode to string

For an Excel Sheet with the structure
Sites   Avail1   Avail2   Avail3
HS1     100      150      175
HS2     200      220      250
..      ...
data from the same row can be kept together in the structure 
{site:[avail1,avail2,avail3], ...}
by using a matrix transpose and calling to_dict with the list argument.
Example 3:
  u_alpha = alpha_df.set_index('Sites').T.to_dict('list')
  alpha = {str(k): v for k, v in u_alpha.items()} # <-- This method can also be used to type-cast values if necessary
  
The more complex the model, the more the input parsing depends on problem formatting and personal preference. 

Right now my method for entering multivariable parameters is to create a dictionary with tuple keys with the format
{(i,j,k):value, ...}
by looping through previously defined sets/imported dictionaries.

For an Excel Sheet with the structure
Sites   PW1   PW2   ..  CS3  CS4
HS1     100   150   ..  175  100
HS2     200   220   ..  250  0.0
..      ...
where feeds and times are appropriately declared sets, the full code to parse it into the form
{(site,feed,time):amount, .. } 

Example 4
  alpha_df = file.parse(header=1, skipfooter=677,parse_cols='B:J')
  u_alpha = alpha_df.set_index('Sites').T.to_dict('list')
  alpha_temp = {str(k): v for k, v in u_alpha.items()}

  alpha = {}
  for k,v in alpha_temp.iteritems():
    i = 0
    for c in feeds:
        for t in times:
            alpha[(k,c,t)] = v[i]
            i = i+1  

