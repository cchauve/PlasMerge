#!/usr/bin/env python

'''
Manipulating pandas dataframes
'''

import pandas as pd

def add_column_from_dict(in_df, in_dict, col_name, reformat_fun=lambda x: x):
    '''
    Add to in_df inplace a column of name col_name from dictionary in_dict 
    with keys being idex of in_df and where each value is transformed by 
    function reformat_fun
    '''
    for idx,_ in in_df.iterrows():
        in_df.at[idx, col_name] = reformat_fun(in_dict[idx])

def add_column_default(in_df, col_name, default_val):
    '''
    Add to in_df inplace a column of name col_name with value default_val
    '''
    for idx,_ in in_df.iterrows():
        in_df.at[idx, col_name] = default_val

def add_column_function(in_df, col_name, col_list, combine_fun):
    '''
    Add to in_df inplace a column of name col_name by applying function
    combine_fun to columns in col_list
    combine_fun is a function that takes as input a list of same length than col_list
    '''
    for idx,row in in_df.iterrows():
        col_vals = [row[col_list[i]] for i in range(len(col_list))]
        in_df.at[idx, col_name] = combine_fun(col_vals)

def modify_column_function(in_df, col_name, modif_fun):
    '''
    Modify in_df[col_name] inplace by applying modif_fun
    '''
    for idx,row in in_df.iterrows():
        prev_val = row[col_name]
        new_val = modif_fun(prev_val)
        in_df.at[idx, col_name] = new_val

