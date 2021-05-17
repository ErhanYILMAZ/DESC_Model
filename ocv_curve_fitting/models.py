"""
Models for the OCV calculations
"""

import pandas as pd


class BatteryScript:
    """
    Script or experiment performed on the battery cell.
    """

    def __init__(self, csvdata):
        """
        Initialize the script measurements.
        """
        columns = ['AhDch', 'AhCha', 'Step', 'VCell1'] 
        #df = pd.read_csv(csvdata, usecols=columns,skiprows=12, sep=';',decimal=',',thousands='.')
        df = pd.read_csv(csvdata)

        #df = df.iloc[1:]
        #df = df.iloc[:-1]


        df['VCell1'] = df['VCell1'].astype(str)
        df['VCell1'] = df['VCell1'].apply(lambda x: x.replace(',', '.')).astype('float')

        df['AhCha'] = df['AhCha'].astype(str)
        df['AhCha'] = df['AhCha'].apply(lambda x: x.replace(',', '.')).astype('float')

        df['AhDch'] = df['AhDch'].astype(str)
        df['AhDch'] = df['AhDch'].apply(lambda x: x.replace(',', '.')).astype('float')

        self.disAh = df['AhDch'].values
        self.chgAh = df['AhCha'].values
        self.step = df['Step'].values
        self.voltage = df['VCell1'].values


class BatteryData:
    """
    Object to store battery measurements from script or experiment for a
    certain temperature.
    """

    def __init__(self, csvfiles):
        """
        Initialize with list of CSV data files.
        """
        self.s1 = BatteryScript(csvfiles[0])
        self.s2 = BatteryScript(csvfiles[1])
        self.s3 = BatteryScript(csvfiles[2])
        self.s4 = BatteryScript(csvfiles[3])


class FileData:
    """
    Calculated data from file.
    """

    def __init__(self, disV, disZ, chgV, chgZ, rawocv, temp):
        self.disV = disV
        self.disZ = disZ
        self.chgV = chgV
        self.chgZ = chgZ
        self.rawocv = rawocv
        self.temp = temp


class ModelOcv:
    """
    Model representing OCV results.
    """

    def __init__(self, OCV0, OCVrel, SOC, OCV, SOC0, SOCrel, OCVeta, OCVQ):
        self.OCV0 = OCV0
        self.OCVrel = OCVrel
        self.SOC = SOC
        self.OCV = OCV
        self.SOC0 = SOC0
        self.SOCrel = SOCrel
        self.OCVeta = OCVeta
        self.OCVQ = OCVQ


