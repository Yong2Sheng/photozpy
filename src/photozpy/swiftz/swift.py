#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 15:00:21 2023

@author: yongsheng
"""

import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord, ICRS, Galactic, FK4, FK5
from swifttools.swift_too import Data, ObsQuery, TOORequests
import time
import shutil
import glob

class UVOTZ():

    def __init__(self, analysis_root_dir, src_catalog_dir = None):
        """Defines the initial parameters.
        Parameters
        ----------
        analysis_root_path: str; the directory where the data will be downloaded and analyzed
        src_catalog_dir: str; the csv file that stores the source name, coordiantes and time window.
        """

        self.analysis_root_dir = analysis_root_dir
        self.full_filter_list = ["ubb", "um2", "uuu", "uvv", "uw1", "uw2"]
        self.src_catalog_dir = src_catalog_dir

        if not self.src_catalog_dir:
            print("You haven't defined the src meta info yet!")
            df = pd.DataFrame(columns=["name", "ra", "dec", "window_lower", "window_upper"])
            df.to_csv(self.analysis_root_dir+"/metadata.csv", sep = ",", index=False)
            print(f"A template csv file is generated at {self.analysis_root_dir}!")

        self.data_dir = self.analysis_root_dir + "/data"
        _ = self.create_folder(self.data_dir)  # here I don't need the returned value

    def str_to_list(self, str):
        return str[1:-2].replace("'", "").split(", ")

    def find_keywords(self, file, *keywords):
        """
        Find the value of the keywords from a fits file. It deals with multi-extension fits files
        so you will get a dictionary of all the requested keywords in each extension.
        "ext_No" and "EXTNAME" are mandatory in order to tag the origin of the keyword values.

        Parameters
        ----------
        file: str; the directory to the fits file
        kewords: str(s); the keywords you want to query

        Example run
        -----------
        find_keywords(file_path, "ASPCORR", "HDUCLAS1")
        Out：
        {'ext_No': [1, 2],
        'EXTNAME': ['bb649482961I', 'bb649505811I'],
        'ASPCORR': ['DIRECT', 'DIRECT'],
        'HDUCLAS1': ['IMAGE', 'IMAGE']}
        """

        with fits.open(file) as hdul:
            ext_nums = len(hdul) - 1

            # initialize the dictionary
            dict_ext = {}
            dict_ext["ext_No"] = [i for i in np.arange(1, ext_nums+1)]
            dict_ext["EXTNAME"] = [None] * ext_nums  # [None, None, None, ...]
            for i in keywords:
                dict_ext[i] = [None] * ext_nums  # generate the keys for the keywords like ext_No and EXTNAME

            
            # read and record the keywords in the dictionary
            n_ = 1
            while n_ <= ext_nums:
                dict_ext["EXTNAME"][n_-1] = hdul[n_].header["EXTNAME"]
                for i in keywords:
                    dict_ext[i][n_-1] = hdul[n_].header[i]
                n_ += 1

        return dict_ext

    def create_folder(self, folder_dir):
        """
        Check if the folder exists and create the folder if not.

        Parameters
        ----------
        folder_path : str; the path of the folder you want to create.

        Returns
        -------
        path_exist : boolean; False if the folder_path doesn't exist;
                     True if the folder_path already exist.
        """

        if os.path.isdir(folder_dir) != True:
            os.mkdir(folder_dir)
            path_exist = False
        else:
            path_exist = True
            
        return path_exist

    def del_then_create_folder(self, folder_dir):
        """
        Delete folder if exists and then create the folder.

        Parameters
        ----------
        folder_dir: str; the directory of the folder
        """
        
        if os.path.exists(folder_dir):
            shutil.rmtree(folder_dir)
        os.mkdir(folder_dir)
        
        return

    def get_obsquery_input(self, **kwargs):
        """
        Generate the input for Swiftttools.
        swifttools accept various types of inputs including obsid, targetid and ra&dec.
        Here the most used one is ra&dec, which is also the default input.

        Parameters
        ----------
        """

        if "obsid" in kwargs.keys():
            return {"name":kwargs["name"],
                    "obsid": kwargs["obsid"]}
            
        elif "targetid" in kwargs.keys():
            return {"name":kwargs["Name"],
                    "targetid": kwargs["targetid"]}

        elif "ra" and "dec" in kwargs.keys():
            ra_dec = {key: value for key, value in kwargs.items() if key in {"ra", "dec"}}
            skycoord = SkyCoord(**ra_dec, unit = (u.hourangle, u.deg), frame = "icrs")
            
            return {"name":kwargs["name"],
                    "skycoord": skycoord}
            
    def download_swift_data(self, radius = 5/60, uvotmode = "0x30ed"):

        # read the csv file
        df = pd.read_csv(self.src_catalog_dir, sep = ",")
        row_nums = df.shape[0]

        obsid = []  # the obsid for all the sources
        obsid_saving_dir = []  # the obsid path for all the sources
        obsid_time = []  # the obsid time for all the sources
        files = []
        filters = []

        for i in np.arange(row_nums):

            info_ = dict(df.iloc[i])  # convert each row into a dict that contains the src info
            src_name = info_["name"].replace(" ", "_")
            src_ra = info_["ra"]
            src_dec = info_["dec"]
            src_window_lower = info_["window_lower"]
            src_window_upper = info_["window_upper"]
            self.src_dir = self.data_dir + f"/{src_name}"
            self.del_then_create_folder(self.src_dir)

            target_info = self.get_obsquery_input(**info_)
            print(target_info)
            oq = ObsQuery(radius = radius, begin = src_window_lower, end = src_window_upper, **target_info)
            id_ = "00000000000"  # this variable will be used to avoid downloading the same data file mutiple times by comparing the observation id

            src_obsid = []
            src_saving_dir = []
            src_obsid_time = []
            for i in np.flip(np.arange(-len(oq), 0)):
                if oq[i].obsid == id_:
                    print(f"{oq[i].obsid} has been downloaded/examined, skipping ...")
                else:
                    if oq[i].uvot_mode == uvotmode:
                        date_ = oq[i].begin.strftime("%Y-%m-%d %H:%M:%S")
                        print(f"{oq[i].obsid} on {date_} is being downloaded, the uvot mode is {oq[i].uvot_mode}")
                        oq[i].download(uvot = True, outdir = self.src_dir, match = ["*/uvot/image/*"])
                        src_obsid.append(oq[i].obsid)  # record the obsid for this source
                        src_saving_dir.append(self.src_dir + "/" + oq[i].obsid)  # record the path to obsids for this source
                        src_obsid_time.append(date_)
                        id_ = oq[i].obsid
                    else:
                        print(f"{oq[i].obsid} uvod mode {oq[i].uvot_mode} is not the one you requested as {uvotmode}, skipping ...")
                        id_ = oq[i].obsid
                    
            obsid.append(src_obsid)
            obsid_saving_dir.append(src_saving_dir)
            obsid_time.append(src_obsid_time)
            print(f"Source {src_name} download finished!")
            print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

        self.meta_data_dir = self.data_dir + "/metadata.csv"
        df["obsid"] = obsid
        df["obsid_saving_dir"] = obsid_saving_dir
        df["obsid_time"] = obsid_time
        df.to_csv(self.meta_data_dir, sep = ",", index=False)

        return

    def check_ASPCORR(self, file):
        """
        Check the ASPCORR keywords.

        Parameters
        ----------
        file: str; the fits file to be exmined.

        Return
        ------
        Boolean
        """
        
        dict_exts = self.find_keywords(file, "ASPCORR")
        asps = dict_exts["ASPCORR"]
        check = ["Yes" for i in asps if "DIRECT" in i]
        if len(check) == len(asps):
            return True
        elif len(check) != len(asps):
            return False

    def get_src_region(self, saving_dir, ra, dec):
        
        reg_files = [saving_dir + "/" + i + ".reg" for i in self.full_filter_list]
        
        for reg in reg_files:
            f = open(reg, "w")
            f.write("# Region file format: DS9 version 4.1\n")
            f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
            f.write("fk5\n")
            f.write(f'circle({ra},{dec},5.000")')
            f.close()
            
        return reg_files

    def sum_images(self, file, sum_type, **kwargs):
        """
        Sum the images from multiple extensions.

        Parameters
        ----------
        file: str; the directory of the file
        sum_type: str; the case of the sum

        Returns
        -------
        out_dir: str; the directory of the output file
        """
        if sum_type == "one_only":  # when there is only one obs file
            filter_name = file[-13:-10]  # get the filer name from the file name
            out_dir = kwargs["summed_dir"] + f"/all{filter_name}.fits"
            os.system(f"uvotimsum infile={file} outfile={out_dir} | tee -a uvotimsum_log.txt >/dev/null 2>&1")

        elif sum_type == "intermediate":  # when there are multiple obs files, sum each by each to get intermediate files
            filter_name = file[-13:-10]
            obsid = file[-24:-13]
            out_dir = kwargs["summed_dir"] + f"/{filter_name}_{obsid}.fits"
            os.system(f"uvotimsum infile={file} outfile={out_dir} | tee -a uvotimsum_log.txt >/dev/null 2>&1")

        elif sum_type == "exclude=NONE":  # when you sum the summed obs files
            filter_name = file[-8:-5]
            out_dir = kwargs["summed_dir"] + f"/all{filter_name}.fits"
            os.system(f"uvotimsum exclude=NONE infile={file} outfile={out_dir} | tee -a uvotimsum_log.txt >/dev/null 2>&1")
            
        else:
            print("Wrong sum type! Please use final or intermediate")
            #break
        return out_dir
        

    def sum_obs_files(self, meta_data_dir = None):
        """
        Sum the multiple obs files for the same source.
        """

        if meta_data_dir == None:
            df = pd.read_csv(self.meta_data_dir, sep = ",")
        else:
            df = pd.read_csv(meta_data_dir, sep = ",")
        df["final_fits"] = None
        df["reg_files"] = None

        for i in np.arange(df.shape[0]):
            final_fits = []
            obsids = self.str_to_list(df.loc[i,"obsid"])
            src_name =  df.loc[i, "name"].replace(" ", "_")
            ra = df.loc[i,"ra"]
            dec = df.loc[i,"dec"]
            print(f"Now summing {df.iloc[i,0]}")

            # create the dir to save the summed images
            summed_dir = self.data_dir + f"/{src_name}" + "/Summed"
            self.create_folder(summed_dir)

            if len(obsids) == 1:  # it means there is only one sky image
                obsid_dir = self.str_to_list(df.loc[i,"obsid_saving_dir"])[0]
                img_files = glob.glob(obsid_dir + "/uvot/image/" + "*sk.img.gz")
                for img_file in img_files:
                    if self.check_ASPCORR(img_file):
                        summed_fits = self.sum_images(img_file, sum_type = "one_only", summed_dir = summed_dir)
                        final_fits.append(summed_fits)

            else:
                for obsid_dir in self.str_to_list(df.loc[i,"obsid_saving_dir"]):
                    img_files = glob.glob(obsid_dir + "/uvot/image/" + "*sk.img.gz")
                    for img_file in img_files:
                        if self.check_ASPCORR(img_file):
                            self.sum_images(img_file, sum_type = "intermediate", summed_dir = summed_dir)
                            
                intermediate_files = glob.glob(summed_dir + "/*fits")  # note fits are the summed sky images
                for filter_ in self.full_filter_list:
                    # find the fits files according to the filer name
                    filter_fits = [fits for fits in intermediate_files if filter_ in fits]
                    if len(filter_fits) == 1:
                        shutil.copy2(filter_fits[0], filter_fits[0][0:-21] + f"/all_{filter_}.fits")

                    else:
                        fappended_file = filter_fits[-1][0:-21] + f"/_{filter_}.fits"  # the file to be appended on
                        shutil.copy2(filter_fits[-1], fappended_file)  # make a copy

                        for j in filter_fits[0:-1]:
                            os.system(f"fappend {j} {fappended_file}")

                        summed_fits = self.sum_images(fappended_file, sum_type="exclude=NONE", summed_dir = summed_dir)
                        final_fits.append(summed_fits)
                        
                _ = self.create_folder(summed_dir + "/intermediate")
                _ = summed_dir + "/*_*"  # only the intermediate files contain _
                inter_dir = summed_dir + "/intermediate"
                os.system(f"mv {_} {inter_dir}")

            reg_files = self.get_src_region(saving_dir = summed_dir, ra = ra, dec = dec)
            
            df.loc[i, "final_fits"] = str(final_fits)
            df.loc[i,"reg_files"] = str(reg_files)

        df.to_csv(self.data_dir+"/metadata.csv", sep = ",", index=False)

    def extract_mag(self, fits_file):
        """
        Extract the AB magnitude from the photometry result file.
        """
        hdul = fits.open(fits_file)
        data = hdul[1].data

        mag = np.round(data["AB_MAG"][0], decimals=2)
        stat_err = np.round(data["AB_MAG_ERR_STAT"][0], decimals=2)
        sys_err = np.round(data["AB_MAG_ERR_SYS"][0], decimals=2)
        err = round(np.sqrt(stat_err**2+sys_err**2),2)
        

        return mag, err

    def extract_filter(self, fits_file):
        """
        Extract the filter info from the photometry result file.
        """
        hdul = fits.open(fits_file)
        FILTER = hdul[1].header["FILTER"]
        if FILTER == "B":
            filter_ = "ubb"
        elif FILTER == "UVM2":
            filter_ = "um2"
        elif FILTER == "U":
            filter_ = "uuu"
        elif FILTER == "V":
            filter_ = "uvv"
        elif FILTER == "UVW1":
            filter_ = "uw1"
        elif FILTER == "UVW2":
            filter_ = "uw2"
        else:
            print("The filter of the image isn't in the filter list!")

        return filter_


    def uvot_photometry(self, meta_data_dir = None):
        """
        Perform photometry.
        """

        if meta_data_dir == None:
            df = pd.read_csv(self.meta_data_dir, sep = ",")
        else:
            df = pd.read_csv(meta_data_dir, sep = ",")

        df_results = pd.DataFrame(columns = ['source_name', "uw2", "uw2_err", "um2", "um2_err", "uw1", "uw1_err", "uuu", "uuu_err", "ubb", "ubb_err", "uvv", "uvv_err"])

        for i in np.arange(df.shape[0]):

            src_name = df.loc[i, "name"].replace(" ", "_")
            
            print(src_name)
            
            dict_new = {'source_name' : [src_name],
                        "uw2" : [-99], 
                        "uw2_err" : [-99],
                        "um2" : [-99],
                        "um2_err" : [-99],
                        "uw1" : [-99],
                        "uw1_err" : [-99],
                        "uuu" : [-99],
                        "uuu_err" : [-99],
                        "ubb" : [-99],
                        "ubb_err" : [-99],
                        "uvv" : [-99],
                        "uvv_err" : [-99]}
            df_new_row = pd.DataFrame.from_dict(dict_new)
            df_results = pd.concat([df_results, df_new_row],ignore_index = True)
            
            os.chdir(f"{self.data_dir}/{src_name}/Summed")
            summed_dir = self.data_dir + f"/{src_name}" + "/Summed"
            
            final_fits = glob.glob("all*.fits")

            for fits_file in final_fits:

                #all_reg_files = self.str_to_list(df.loc[0,"reg_files"])

                # defining all kinds of inputs and outputs also the column names
                filter_ = self.extract_filter(fits_file)
                filter_err = filter_ + "_err"
                src_region_file = summed_dir + f"/{filter_}.reg"
                print(src_region_file)
                bkg_region_file = summed_dir + f"/bg{filter_}.reg"
                outfits = summed_dir + f"/{filter_}_Results.fits"
                outtxt = summed_dir + f"/{filter_}_Results.txt"

                # run the command
                print(f"Running uvotsource image={fits_file} srcreg={src_region_file} bkgreg={bkg_region_file} sigma=3 cleanup=y clobber=y outfile={outfits} | tee {outtxt} >/dev/null")
                os.system(f"uvotsource image={fits_file} srcreg={src_region_file} bkgreg={bkg_region_file} sigma=3 cleanup=y clobber=y outfile={outfits} | tee {outtxt} >/dev/null")

                # extract the magnitudes
                ab_mag, ab_err = self.extract_mag(outfits)

                print(filter_)
                print(ab_mag)
                print(ab_err)
                print("++++++++++++++++++++++++++++++++++++++")

                index = df_results.shape[0]-1  # always append the data to the last row
                df_results.loc[index, filter_] = ab_mag
                df_results.loc[index, filter_err] = ab_err

                df_results.to_csv(self.data_dir+"/Magnitudes.csv", sep = ",", index=False)
                
            os.chdir(self.data_dir)
        
                
                
            
            
        
        
        
                
        
    
            
                        
                        
        
        
        
    
