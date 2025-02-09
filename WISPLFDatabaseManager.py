import datetime
import numpy as np
import sqlite3


class WISPLFDatabaseManager:
    """
    Utility class to manage persistence of line identification data during the
    WISP line-finding process.
    """

    validMutableFlags = [
        'ZEROTH',
        'CONTIN',
        'MISC'
    ]

    validFlags = [
        'CONTAM',
        'REJECT'
    ] + validMutableFlags

    # validContamFlags = ['lya',
    #                     'c4',
    #                     'he2',
    #                     'o3s',
    #                     'c3',
    #                     'mg2',
    #                     'o2',
    #                     'hg',
    #                     'hb',
    #                     'o3',
    #                     'o1',
    #                     'ha',
    #                     's2',
    #                     's31',
    #                     's32',
    #                     'he1',
    #                     'c'] # M.D.R 01/06/2021

    validContamFlags = \
    ['la_1216',
    'n5_1238',
    'n5_1242',
    'c4_1548',
    'c4_1550',
    'h2_1640',
    'o3_1660',
    'o3_1666',
    's3_1883',
    's3_1892',
    'c3_1907',
    'c3_1909',
    'm2_2796',
    'm2_2803',
    'o2_3727',
    'o2_3730',
    'hg_4342',
    'o3_4363',
    'h2_4686',
    'hb_4863',
    'o3_4959',
    'o3_5007',
    'o1_6300',
    'o1_6363',
    'n2_6550',
    'ha_6565',
    'n2_6585',
    's2_6716',
    's2_6731',
    's3_9069',
    's3_9532',
    'he10830',
    'cont'] # MDR 2022/07/22

    tableNames = ['catalogue',
                  'annotations',
                  'flags'
                  ]

    # lineListHeadings = [
    #     'ParNo',
    #     'ObjID',
    #     'RA',
    #     'Dec',
    #     'Jmagnitude [99.0 denotes no detection]',
    #     'Hmagnitude [99.0 denotes no detection]',
    #     'A_IMAGE',
    #     'B_IMAGE',
    #     'redshift',
    #     'redshift_err',
    #     'snr_tot_others', # MDR 2022/05/25
    #     'dz_oiii',
    #     'dz_oii',
    #     'dz_siii_he1',
    #     'G141_FWHM_Obs [Angs]',
    #     'G141_FWHM_Obs_err',
    #     'lya_1216_flux', # M.D.R 01/06/2021
    #     'lya_1216_err', # M.D.R 01/06/2021
    #     'lya_1216_EW_obs', # M.D.R 01/06/2021
    #     'lya_1216_contam', # M.D.R 01/06/2021
    #     'civ_flux',
    #     'civ_error',
    #     'civ_EW_obs',
    #     'civ_contam',
    #     'heii_flux',
    #     'heii_error',
    #     'heii_EW_obs',
    #     'heii_contam',
    #     'oiii_1663_flux',
    #     'oiii_1663_error',
    #     'oiii_1663_EW_obs',
    #     'oiii_1663_contam',
    #     'ciii_flux',
    #     'ciii_error',
    #     'ciii_EW_obs',
    #     'ciii_contam',
    #     'mgii_flux',
    #     'mgii_error',
    #     'mgii_EW_obs',
    #     'mgii_contam',
    #     'oii_flux',
    #     'oii_error',
    #     'oii_EW_obs',
    #     'oii_contam',
    #     'hg_flux',
    #     'hg_err',
    #     'hg_EW_obs',
    #     'hg_contam',
    #     'hb_flux',
    #     'hb_err',
    #     'hb_EW_obs',
    #     'hb_contam',
    #     'oiii_flux [both lines]',
    #     'oiii_err [both lines]',
    #     'oiii_EW_obs [both lines]',
    #     'oiii_contam [both lines]',
    #     'oi_flux',
    #     'oi_error',
    #     'oi_EW_obs',
    #     'oi_contam',
    #     'hanii_flux',
    #     'hanii_err',
    #     'hanii_EW_obs',
    #     'hanii_contam',
    #     'sii_flux',
    #     'sii_err',
    #     'sii_EW_obs',
    #     'sii_contam',
    #     'siii_9069_flux',
    #     'siii_9069_err',
    #     'siii_9069_EW_obs',
    #     'siii_9069_contam',
    #     'siii_9532_flux',
    #     'siii_9532_err',
    #     'siii_9532_EW_obs',
    #     'siii_9532_contam',
    #     'he1_10830_flux',
    #     'he1_10830_err',
    #     'he1_10830_EW_obs',
    #     'he1_10830_contam',
    # ]

    # fitResultKeys = ['redshift',
    #                  'redshift_err',
    #                  'dz_oiii',
    #                  'dz_oii',
    #                  'dz_siii_he1',
    #                  'fwhm_g141',
    #                  'fwhm_g141_err',
    #                  'oii_flux',
    #                  'oii_error',
    #                  'oii_ew_obs',
    #                  'hg_flux',
    #                  'hg_error',
    #                  'hg_ew_obs',
    #                  'hb_flux',
    #                  'hb_error',
    #                  'hb_ew_obs',
    #                  'oiii_flux',
    #                  'oiii_error',
    #                  'oiii_ew_obs',
    #                  'hanii_flux',
    #                  'hanii_error',
    #                  'hanii_ew_obs',
    #                  'sii_flux',
    #                  'sii_error',
    #                  'sii_ew_obs',
    #                  'siii_9069_flux',
    #                  'siii_9069_error',
    #                  'siii_9069_ew_obs',
    #                  'siii_9532_flux',
    #                  'siii_9532_error',
    #                  'siii_9532_ew_obs',
    #                  'he1_flux',
    #                  'he1_error',
    #                  'he1_ew_obs',
    #                  'lya_flux', # M.D.R 01/06/2021
    #                  'lya_error', # M.D.R 01/06/2021
    #                  'lya_ew_obs'] # M.D.R 01/06/2021

    def __init__(self, dbFileNamePrefix):
        #        print('Using sqlite3 version {}'.format(sqlite3.version))
        self.dbFileNamePrefix = dbFileNamePrefix
        self.dbFileName = '{}_sqlite.db'.format(self.dbFileNamePrefix)
        self.dbConnection = sqlite3.connect(self.dbFileName)
        self.dbConnection.row_factory = sqlite3.Row
        self.dbCursor = self.dbConnection.cursor()
        self.checkAndInitTables()

    def __del__(self):
        self.dbConnection.commit()
        self.dbConnection.close()

    def checkAndInitTables(self):
        self.createCatalogueTable()
        self.createAnnotationTable()
        self.createFlagTable()

    def createCatalogueTable(self):
        #        print('Creating catalogue table in SQLLite database...')
        self.dbCursor.execute('''CREATE TABLE IF NOT EXISTS catalogue (
                              ParNo int,
                              ObjID int,
                              RA real,
                              Dec real,
                              Jmagnitude real,
                              Hmagnitude real,
                              A_IMAGE real,
                              B_IMAGE real,
                              redshift real,
                              redshift_err real,
                              dz_oiii real,
                              dz_oii real,
                              dz_siii_he1 real,
                              fwhm_g141 real,
                              fwhm_g141_err real,
                              oii_flux real,
                              oii_error real,
                              oii_ew_obs real,
                              hg_flux real,
                              hg_error real,
                              hg_ew_obs real,
                              hb_flux real,
                              hb_error real,
                              hb_ew_obs real,
                              oiii_flux real,
                              oiii_error real,
                              oiii_ew_obs real,
                              hanii_flux real,
                              hanii_error real,
                              hanii_ew_obs real,
                              sii_flux real,
                              sii_error real,
                              sii_ew_obs real,
                              siii_9069_flux real,
                              siii_9069_error real,
                              siii_9069_ew_obs real,
                              siii_9532_flux real,
                              siii_9532_error real,
                              siii_9532_ew_obs real,
                              he1_flux real,
                              he1_error real,
                              he1_ew_obs real,
                              lya_flux real,
                              lya_error real,
                              lya_ew_obs real,
                              ContamFlag int,
                              EntryTime text
                              )''')
        self.dbConnection.commit() # M.D.R 01/06/2021 - added lya above
#        print('Done.' if self.dbCursor.rowcount >
#              0 else 'Table was already present.')

    def createAnnotationTable(self):
        #        print('Creating annotations table in SQLLite database...')
        self.dbCursor.execute('''CREATE TABLE IF NOT EXISTS annotations (
                         ParNo int,
                         ObjID int,
                         Comment text
                         )''')
        self.dbConnection.commit()
#        print('Done.' if self.dbCursor.rowcount >
#              0 else 'Table was already present.')

    def createFlagTable(self):
        #        print('Creating annotations table in SQLLite database...')
        self.dbCursor.execute('''CREATE TABLE IF NOT EXISTS flags (
                         ParNo int,
                         ObjID int,
                         FlagName text,
                         FlagValue int
                         )''')
        self.dbConnection.commit()
#        print('Done.' if self.dbCursor.rowcount >
#              0 else 'Table was already present.')

    def saveCatalogueEntry(self, catalogueEntryData):
        query = 'INSERT INTO catalogue VALUES ({})'.format(
            ','.join(['?'] * len(catalogueEntryData))
        )
        # print(query)
        self.dbCursor.execute(query, catalogueEntryData)
        self.dbConnection.commit()

    def loadCatalogueEntry(self, parNumber, objectId):
        query = 'SELECT * FROM catalogue WHERE (ParNo=? AND ObjID=?)'
        self.dbCursor.execute(query, (parNumber, objectId))
        catalogueEntryData = self.dbCursor.fetchall()
        if len(catalogueEntryData) < 1:
            return None
        # FIXME: For now, just return the last row added
        nonFitResults = [catalogueEntryData[-1][key]
                         for key in list(catalogueEntryData[-1].keys()) if key not in WISPLFDatabaseManager.fitResultKeys]
        fitResults = {key: catalogueEntryData[-1][key]
                      for key in WISPLFDatabaseManager.fitResultKeys}
        return tuple(nonFitResults), fitResults

    def getMostRecentObject(self, parNumber=None):
        query = '''SELECT ParNo, ObjID, EntryTime
        FROM catalogue{}
        ORDER BY DATETIME(EntryTime)
        DESC LIMIT 1'''.format('WHERE ParNo = ?' if parNumber is not None else '')
        self.dbCursor.execute(
            *[argument for argument in [query, parNumber] if argument is not None])
        mostRecentEntry = self.dbCursor.fetchone()
        if mostRecentEntry is not None:
            return mostRecentEntry['ParNo'], mostRecentEntry['ObjId'], mostRecentEntry['EntryTime']
        return None

    def layoutCatalogueData(self,
                            parNumber,
                            objectId,
                            ra,
                            dec,
                            jMagnitude,
                            hMagnitude,
                            aImage,
                            bImage,
                            fitResults,
                            flagContent):
        if fitResults is None:
            fitResults = {
                key: None for key in WISPLFDatabaseManager.fitResultKeys}
        return (parNumber, objectId, ra, dec, jMagnitude, hMagnitude, aImage, bImage,
                fitResults['redshift'],
                fitResults['redshift_err'],
                fitResults['dz_oiii'],
                fitResults['dz_oii'],
                fitResults['dz_siii_he1'],
                fitResults['fwhm_g141'],
                fitResults['fwhm_g141_err'],
                fitResults['oii_flux'],
                fitResults['oii_error'],
                fitResults['oii_ew_obs'],
                fitResults['hg_flux'],
                fitResults['hg_error'],
                fitResults['hg_ew_obs'],
                fitResults['hb_flux'],
                fitResults['hb_error'],
                fitResults['hb_ew_obs'],
                fitResults['oiii_flux'],
                fitResults['oiii_error'],
                fitResults['oiii_ew_obs'],
                fitResults['hanii_flux'],
                fitResults['hanii_error'],
                fitResults['hanii_ew_obs'],
                fitResults['sii_flux'],
                fitResults['sii_error'],
                fitResults['sii_ew_obs'],
                fitResults['siii_9069_flux'],
                fitResults['siii_9069_error'],
                fitResults['siii_9069_ew_obs'],
                fitResults['siii_9532_flux'],
                fitResults['siii_9532_error'],
                fitResults['siii_9532_ew_obs'],
                fitResults['he1_flux'],
                fitResults['he1_error'],
                fitResults['he1_ew_obs'],
                fitResults['lya_flux'], # M.D.R 01/06/2021
                fitResults['lya_error'], # M.D.R 01/06/2021
                fitResults['lya_ew_obs'], # M.D.R 01/06/2021
                flagContent,
                str(datetime.datetime.now().isoformat())
                )

    def saveAnnotation(self, annotationData):
        query = 'INSERT INTO annotations VALUES ({})'.format(
            ','.join(['?'] * len(annotationData))
        )
        self.dbCursor.execute(query, annotationData)
        self.dbConnection.commit()

    def setFlagsFromString(self, parNumber, objId, flagDataString, delimiter=','):
        # Assumes that string is a delimiter separated list of flags and values
        flagDataTokens = [token.strip()
                          for token in flagDataString.split(delimiter)]
        flagNames = flagDataTokens[::2]
        flagValues = flagDataTokens[1::2]
        flagData = list(zip(flagNames, flagValues))
        self.setFlags(parNumber, objId, flagData)

    def setFlags(self, parNumber, objId, flagData):
        # Only attempt to set values for valid flags
        flagData = [(parNumber, objId, flagDatum[0], flagDatum[1])
                    for flagDatum in flagData
                    if flagDatum[0] in WISPLFDatabaseManager.validFlags + WISPLFDatabaseManager.validContamFlags ]
        # flagData can be a list of (flagName, flagValue) tuples
        query = 'INSERT INTO flags VALUES (?, ?, ?, ?)'
        self.dbCursor.executemany(query, flagData)
        self.dbConnection.commit()

    def resetDatabaseTables(self):
        query = 'DROP TABLE IF EXISTS {}'
        for tableName in WISPLFDatabaseManager.tableNames :
            self.dbCursor.execute(query.format(tableName))
            self.dbConnection.commit()
        self.checkAndInitTables()

    def writeCatalogueTextFile(self):
        catalogueQuery = 'SELECT * FROM catalogue'
