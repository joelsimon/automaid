import glob
import os
import shutil
import json
import requests
import struct
import re
import traceback

mermaid_path = os.environ["MERMAID"]
database_path = os.path.join(mermaid_path, "database")

DATABASE_LINK_NAME = "Databases.json"
def database_update(path):
    print("Update Databases")
    network = 1
    database_path = path
    try:
        # Get linker file (link database and profiler version)
        url = 'http://mermaid.osean.fr/databases/' + DATABASE_LINK_NAME
        request = requests.get(url, auth=('osean', 'osean3324'), timeout=10)
    except Exception as e:
        print("Exception: \"" + str(e) + "\" detected when get " + DATABASE_LINK_NAME)
        network = 0
    else:
        if request.status_code == 200:
            # Read link file content into object
            database_list = request.json()
            # Retreive list of databases into linker file
            for database in database_list:
                if database["Name"]:
                    try:
                        # Get Database files (one by one)
                        url_database = 'http://mermaid.osean.fr/databases/' + database["Name"]
                        new_req = requests.get(url_database, auth=('osean', 'osean3324'), timeout=10)
                        # Store temporary database into linker object
                        database["data"] = new_req.json()
                        if new_req.status_code != 200:
                            print("Error " + str(new_req.status_code) + " when get " + database["Name"])
                            network = 0
                    except Exception as e:
                        print("Exception: \"" + str(e) + "\" detected when get " + str(database["Name"]))
                        network = 0
        else:
            print("Error " + str(request.status_code) + " when get " + DATABASE_LINK_NAME)
            network = 0
        if network > 0:
            # All databases received correctly => delete current files
            if os.path.exists(database_path):
                shutil.rmtree(database_path)
            os.makedirs(database_path)
            # Store databases content
            for database in database_list:
                if database["Name"]:
                    path = os.path.join(database_path, database["Name"])
                    with open(path, 'w') as databasefile:
                        json.dump(database["data"], databasefile)
                    # Delete content from link object
                    database.pop('data', None)
            # Store link file
            link_path = os.path.join(database_path, DATABASE_LINK_NAME)
            with open(link_path, 'w') as linkfile:
                json.dump(database_list, linkfile, indent=4)

# Get database name with linker file and version read on file
def database_get_version(file_version,model) :
    link_path = os.path.join(database_path,DATABASE_LINK_NAME)
    if os.path.exists(link_path):
        with open(link_path,"r") as f:
            databases = json.loads(f.read())
        # get major and minor versions
        file_version=file_version.split(".")
        file_major = 2
        file_minor = 17
        if file_version[0] :
            file_major = int(file_version[0])
        if file_version[1] :
            file_minor = int(file_version[1])
        for database in databases :
            if not database["Model"] or (model == database["Model"]):
                database_minor_max = 2147483647
                database_major_max = 2147483647
                databaseMin_version = database["MinVersion"].split(".")
                database_major_min = int(databaseMin_version[0])
                database_minor_min = int(databaseMin_version[1])
                if database["MaxVersion"] != "None":
                    databaseMax_version = database["MaxVersion"].split(".")
                    database_major_max = int(databaseMax_version[0])
                    database_minor_max = int(databaseMax_version[1])
                if ((file_major > database_major_min) and (file_major < database_minor_max)) \
                or ((file_major >= database_major_min) and (file_major <= database_major_max) \
                and (file_minor >= database_minor_min) and (file_minor <= database_minor_max)) :
                    return database["Name"]
                #print "\r\n"
        print(("No Database available for : " + str(file_version)))
        return ""
    else :
        print(("No link file : " + link_path))
        return ""

 # Concatenate .[0-9][0-9][0-9] files .LOG and .BIN files in the path
def concatenate_files(path):
    # List log files (not cyphered)
    files_path = glob.glob(path + "*.LOG")
    # Add bin files
    files_path += glob.glob(path + "*.BIN")
    # Concatenate all files
    for file_path in files_path:
        bin = b''
        # list file with same head than log file
        files_to_merge = list(glob.glob(file_path[:-4] +".[0-9][0-9][0-9]"))
        # Add end file to the list
        files_to_merge.append(file_path)
        # Sort list => [0000_XXXXXX.000,0000_XXXXXX.001,0000_XXXXXX.002,0000_XXXXXX.LOG]
        files_to_merge.sort()
        if len(files_to_merge) > 1:
            # Need to merge multiple file
            print("Files to merge :" + repr(files_to_merge))
            for file_to_merge in files_to_merge :
                if file_to_merge[-3:].isdigit():
                    # Extension with digit => append to string
                    with open(file_to_merge, "rb") as fl:
                        bin += fl.read()
                    # Remove file after read
                    os.remove(file_to_merge)
                else :
                    if len(bin) > 0:
                        # If log extension is not a digit and the log string is not empty
                        # we need to add it at the end of the file
                        with open(file_to_merge, "rb") as fl:
                            bin += fl.read()
                        with open(file_to_merge, "wb") as fl:
                            fl.write(bin)
                        bin = b''

def decrypt_explicit(f,LOG_card,WARN_card,ERR_card) :
    #Read head
    string = ""
    IDbytes = f.read(2)
    if len(IDbytes) != 2 :
        print("err:IDbytes")
        return "err:IDbytes\r\n";
    TIMESTAMPbytes = f.read(4)
    if len(TIMESTAMPbytes) != 4 :
        print("err:TIMESTAMPbytes")
        return "err:TIMESTAMPbytes\r\n";
    INFOSbytes = f.read(1)
    if INFOSbytes == "" :
        print("err:INFOSbytes")
        return "err:INFOSbytes\r\n";
    DATASIZEbytes = f.read(1)
    if DATASIZEbytes == "" :
        print("err:DATASIZEbytes")
        return "err:DATASIZEbytes\r\n";
    #unpack head
    try :
        id = struct.unpack('<H', IDbytes)[0]
        timestamp = struct.unpack('<I', TIMESTAMPbytes)[0]
        infos = struct.unpack('<B', INFOSbytes)[0]
        dataSize = struct.unpack('<B', DATASIZEbytes)[0]
    except :
        traceback.print_exc()
        print("err:header")
        return "err:header\r\n";

    #Process head
    idString = "0x"+"{0:0{1}X}".format(id,4)+"UL"
    binaryinfo = "{0:08b}".format(infos)
    logtype = "00"
    argformat = "00"
    logtype = binaryinfo[-2:]
    argformat = binaryinfo[-4:-2]
    if argformat != "00":
        return "err:argformat\r\n";

    #print ("ID : " + str(id))
    #print ("IDString : " + str(idString))
    #print ("Timestamp : " + str(timestamp))
    #print ("Infos : " + str(infos))
    #print ("BinaryInfos : " + str(binaryinfo))
    #print ("Type : " + str(type))
    #print ("ArgFormat : " + str(argformat))
    #print ("dataSize : " + str(dataSize))

    decrypt_card={}
    type_string = ""
    if logtype == "00":
        decrypt_card = LOG_card
    elif logtype == "01":
        type_string = "<WARN>"
        decrypt_card = WARN_card
    elif logtype == "10":
        type_string = "<ERR>"
        decrypt_card = ERR_card
    else :
        type_string = "<DBG>"

    Formats = []
    File = "MAIN"
    Level = "1"

    if(id < len(decrypt_card)) :
        index = id
    else :
        index = len(decrypt_card)-1
    while index >= 0 :
        if decrypt_card[index]["ID"] == idString:
            Formats = decrypt_card[index]["FORMATS"]
            File = decrypt_card[index]["FILE"]
            Level = decrypt_card[index]["LEVEL"]
            break;
        index = index - 1

    if len(Formats) <= 0 :
        f.read(dataSize)
        string += str(timestamp) + ":" + type_string + "["+"{:04d}".format(id)+"] Format not found\r\n"
        return string
    #print (Formats)
    string += str(timestamp) + ":" + type_string
    string +="["+"{:6}".format(File)+","+"{:04d}".format(id)+"]"
    index=0
    argIndex=0
    if dataSize > 0:
        while index < dataSize :
            #Read Argument Head
            ARGINFOSByte = f.read(1)
            if ARGINFOSByte == "":
                print("err:ARGINFOSByte")
                return "err:ARGINFOSByte\r\n"
            ARGSIZEByte = f.read(1)
            if ARGSIZEByte == "":
                print("err:ARGSIZEByte")
                return "err:ARGSIZEByte\r\n"
            #Unpack Argument Head
            try :
                ArgInfos=struct.unpack('<B', ARGINFOSByte)[0]
                ArgSize = struct.unpack('<B', ARGSIZEByte)[0]
            except :
                traceback.print_exc()
                print("err:INFOSSIZE")
                return "err:INFOSSIZE\r\n"

            #Process Argument Head
            ArgInfosBinary="{0:08b}".format(ArgInfos)
            ArgType = ArgInfosBinary[-2:]

            #print ("ArgInfosBinary : " + str(ArgInfosBinary))
            #print ("ArgType : " + str(ArgType))
            #print ("ArgSize : " + str(ArgSize))
            index = index+2
            Formats[argIndex] = Formats[argIndex].replace(r"\r\n","\r\n")
            if ArgSize > 0:
                Arg = 0
                if ArgType == "00":
                    if ArgSize == 4:
                        ArgByte = f.read(4)
                        if len(ArgByte) != 4 :
                            print("err:TYPE00SIZE04")
                            return "err:TYPE00SIZE04\r\n"
                        try :
                            Arg = struct.unpack('<i', ArgByte)[0]
                        except :
                            traceback.print_exc()
                            print("err:TYPE00SIZE04")
                            return "err:TYPE00SIZE04\r\n"
                    elif ArgSize == 2:
                        ArgByte = f.read(2)
                        if len(ArgByte) != 2 :
                            print("err:TYPE00SIZE02")
                            return "err:TYPE00SIZE02\r\n"
                        try :
                            Arg = struct.unpack('<h', ArgByte)[0]
                        except :
                            traceback.print_exc()
                            print("err:TYPE00SIZE02")
                            return "err:TYPE00SIZE02\r\n"
                    elif ArgSize == 1:
                        ArgByte = f.read(1)
                        if ArgByte == "" :
                            print("err:TYPE00SIZE01")
                            return "err:TYPE00SIZE01\r\n"
                        try :
                            Arg = struct.unpack('<b', ArgByte)[0]
                        except :
                            traceback.print_exc()
                            print("err:TYPE00SIZE01")
                            return "err:TYPE00SIZE01\r\n"
                elif ArgType == "01":
                    # unsigned integer
                    if ArgSize == 4:
                        ArgByte = f.read(4)
                        if len(ArgByte) != 4 :
                            print("err:TYPE01SIZE04")
                            return "err:TYPE01SIZE04\r\n"
                        try :
                            Arg = struct.unpack('<I', ArgByte)[0]
                        except :
                            traceback.print_exc()
                            print("err:TYPE01SIZE04")
                            return "err:TYPE01SIZE04\r\n"
                    elif ArgSize == 2:
                        ArgByte = f.read(2)
                        if len(ArgByte) != 2 :
                            print("err:TYPE01SIZE02")
                            return "err:TYPE01SIZE02\r\n"
                        try :
                            Arg = struct.unpack('<H', ArgByte)[0]
                        except :
                            traceback.print_exc()
                            print("err:TYPE01SIZE02")
                            return "err:TYPE01SIZE02\r\n"
                    elif ArgSize == 1:
                        ArgByte = f.read(1)
                        if ArgByte == "":
                            print("err:TYPE01SIZE01")
                            return "err:TYPE01SIZE01\r\n"
                        try :
                            Arg = struct.unpack('<B', ArgByte)[0]
                        except :
                            traceback.print_exc()
                            print("err:TYPE01SIZE01")
                            return "err:TYPE01SIZE01\r\n"
                elif ArgType == "11":
                    # string
                    ArgByte = f.read(ArgSize)
                    if len(ArgByte) != ArgSize :
                        print("err:TYPE11")
                        return "err:TYPE11\r\n"
                    if ArgByte[ArgSize-1] == 0 :
                        ArgByte = ArgByte[:-1]
                    Arg = ArgByte
                    try :
                        Arg = Arg.decode('ascii', 'ignore')
                    except :
                        traceback.print_exc()
                        print("err:TYPE11")
                        return "err:TYPE11\r\n"
                    #replace none ascii characters
                #print (ArgSize)
                #print ("Format : " + str(Formats[argIndex]) + "\r\n")
                try :
                    if "%.*s" in Formats[argIndex]:
                        string += (Formats[argIndex] % (ArgSize,Arg))
                    else :
                        string += (Formats[argIndex] % Arg)
                except :
                    traceback.print_exc()
                    print("error format")
                    return "err:format\r\n"
            else :
                #print ("Format : " + str(Formats[argIndex]) + "\r\n")
                string += str(Formats[argIndex])
                index = index + 1
            index = index + ArgSize
            argIndex = argIndex + 1
    else :
        #print ("Format : " + (Formats[0].replace(r"\r\n","\r\n")) + "\r\n")
        string += str(Formats[0].replace(r"\r\n","\r\n"))
    string += "\r\n"
    return string


def decrypt_short(f) :
    string = ""
    id = f.read(1)
    format = ""
    shortId = 256
    if id == "" :
        print("err:EMPTYSHORTID")
        return "err:EMPTYSHORTID\r\n";
    try :
        shortId = struct.unpack('<B', id)[0]
    except :
        traceback.print_exc()
        print("err:UNPACKSHORTID")
        return "err:UNPACKSHORTID\r\n"
    if shortId == 0 :
        # Pressure LOG
        TIMESTAMPbytes = f.read(4)
        if len(TIMESTAMPbytes) != 4 :
            print("err:TIMESTAMPbytes")
            return "err:TIMESTAMPbytes\r\n";
        # pressure measure
        PRESSbytes = f.read(4)
        if len(PRESSbytes) != 4 :
            print("err:PRESSbytes")
            return "err:PRESSbytes\r\n";
        timestamp = 0
        pressure = 0
        try :
            # Unpack Integer of 4 bytes
            timestamp = struct.unpack('<I', TIMESTAMPbytes)[0]
            pressure = struct.unpack('<l', PRESSbytes)[0]
        except :
            traceback.print_exc()
            print("err:UNPACKPRESS")
            return "err:UNPACKPRESS\r\n"
        format = str(timestamp) + ":[PRESS ,0038]P%+7dmbar\r\n"
        return format % pressure
    elif shortId == 1 :
        # Pump time LOG
        TIMESTAMPbytes = f.read(4)
        if len(TIMESTAMPbytes) != 4 :
            print("err:TIMESTAMPbytes")
            return "err:TIMESTAMPbytes\r\n";
        PUMPTIMEbytes = f.read(4)
        if len(PUMPTIMEbytes) != 4 :
            print("err:PUMPTIMEbytes")
            return "err:PUMPTIMEbytes\r\n";
        timestamp = 0
        pump_time = 0
        try :
            # Unpack Integer of 4 bytes
            timestamp = struct.unpack('<I', TIMESTAMPbytes)[0]
            pump_time = struct.unpack('<l', PUMPTIMEbytes)[0]
        except :
            traceback.print_exc()
            print("err:UNPACKPUMPTIME")
            return "err:UNPACKPUMPTIME\r\n"
        format = str(timestamp) + ":[PUMP  ,0016]pump during %dms\r\n"
        return format % pump_time
    elif shortId == 2 :
        # Valve time LOG
        TIMESTAMPbytes = f.read(4)
        if len(TIMESTAMPbytes) != 4 :
            print("err:TIMESTAMPbytes")
            return "err:TIMESTAMPbytes\r\n";
        VALVETIMEbytes = f.read(4)
        if len(VALVETIMEbytes) != 4 :
            print("err:VALVETIMEbytes")
            return "err:VALVETIMEbytes\r\n";
        timestamp = 0
        valve_time = 0
        try :
            # Unpack Integer of 4 bytes
            timestamp = struct.unpack('<I', TIMESTAMPbytes)[0]
            valve_time = struct.unpack('<l', VALVETIMEbytes)[0]
        except :
            traceback.print_exc()
            print("err:UNPACKVALVETIME")
            return "err:UNPACKVALVETIME\r\n"
        format = str(timestamp) + ":[VALVE ,0034]valve opening %dms\r\n"
        return format % valve_time
    elif shortId == 3 :
        # Bypass time LOG
        TIMESTAMPbytes = f.read(4)
        if len(TIMESTAMPbytes) != 4 :
            print("err:TIMESTAMPbytes")
            return "err:TIMESTAMPbytes\r\n";
        BYPASSTIMEbytes = f.read(4)
        if len(BYPASSTIMEbytes) != 4 :
            print("err:BYPASSTIMEbytes")
            return "err:BYPASSTIMEbytes\r\n";
        timestamp = 0
        bypass_time = 0
        try :
            # Unpack Integer of 4 bytes
            timestamp = struct.unpack('<I', TIMESTAMPbytes)[0]
            bypass_time = struct.unpack('<l', BYPASSTIMEbytes)[0]
        except :
            traceback.print_exc()
            print("err:UNPACKBYPASSTIME")
            return "err:UNPACKBYPASSTIME\r\n"
        format = str(timestamp) + ":[BYPASS,0035]bypass opening %dms\r\n"
        return format % bypass_time
    elif shortId == 4 :
        # Pump power LOG
        TIMESTAMPbytes = f.read(4)
        if len(TIMESTAMPbytes) != 4 :
            print("err:TIMESTAMPbytes")
            return "err:TIMESTAMPbytes\r\n";
        PUMPVOLTbytes = f.read(4)
        if len(PUMPVOLTbytes) != 4 :
            print("err:PUMPVOLTbytes")
            return "err:PUMPVOLTbytes\r\n";
        PUMPCURRENTbytes = f.read(4)
        if len(PUMPCURRENTbytes) != 4 :
            print("err:PUMPCURRENTbytes")
            return "err:PUMPCURRENTbytes\r\n";
        PUMPPRESSbytes = f.read(4)
        if len(PUMPPRESSbytes) != 4 :
            print("err:PUMPPRESSbytes")
            return "err:PUMPPRESSbytes\r\n";
        timestamp = 0
        pump_volt = 0
        pump_current = 0
        pump_press = 0
        try :
            # Unpack Integer of 4 bytes
            timestamp = struct.unpack('<I', TIMESTAMPbytes)[0]
            pump_volt = struct.unpack('<l', PUMPVOLTbytes)[0]
            pump_current = struct.unpack('<l', PUMPCURRENTbytes)[0]
            pump_press = struct.unpack('<l', PUMPPRESSbytes)[0]
        except :
            traceback.print_exc()
            print("err:UNPACKPUMPPOWER")
            return "err:UNPACKPUMPPOWER\r\n"
        format = str(timestamp) + ":[PUMP  ,0212]battery %5dmV, %7duA (steady state), P%+7dmbar\r\n"
        return format % (pump_volt,pump_current,pump_press)
    elif shortId == 5 :
        # Valve power LOG
        TIMESTAMPbytes = f.read(4)
        if len(TIMESTAMPbytes) != 4 :
            print("err:TIMESTAMPbytes")
            return "err:TIMESTAMPbytes\r\n";
        VALVEVOLTbytes = f.read(4)
        if len(VALVEVOLTbytes) != 4 :
            print("err:VALVEVOLTbytes")
            return "err:VALVEVOLTbytes\r\n";
        VALVECURRENTbytes = f.read(4)
        if len(VALVECURRENTbytes) != 4 :
            print("err:VALVECURRENTbytes")
            return "err:VALVECURRENTbytes\r\n";
        VALVEPRESSbytes = f.read(4)
        if len(VALVEPRESSbytes) != 4 :
            print("err:VALVEPRESSbytes")
            return "err:VALVEPRESSbytes\r\n";
        timestamp = 0
        valve_volt = 0
        valve_current = 0
        valve_press = 0
        try :
            # Unpack Integer of 4 bytes
            timestamp = struct.unpack('<I', TIMESTAMPbytes)[0]
            valve_volt = struct.unpack('<l', VALVEVOLTbytes)[0]
            valve_current = struct.unpack('<l', VALVECURRENTbytes)[0]
            valve_press = struct.unpack('<l', VALVEPRESSbytes)[0]
        except :
            traceback.print_exc()
            print("err:UNPACKVALVEPOWER")
            return "err:UNPACKVALVEPOWER\r\n"
        format = str(timestamp) + ":[VALVE ,0234]battery %5dmV, %7duA, P%+7dmbar\r\n"
        return format % (valve_volt,valve_current,valve_press)
    elif shortId == 6 :
        # Bpopen power LOG
        TIMESTAMPbytes = f.read(4)
        if len(TIMESTAMPbytes) != 4 :
            print("err:TIMESTAMPbytes")
            return "err:TIMESTAMPbytes\r\n";
        BPOPENVOLTbytes = f.read(4)
        if len(BPOPENVOLTbytes) != 4 :
            print("err:BPOPENVOLTbytes")
            return "err:BPOPENVOLTbytes\r\n";
        BPOPENCURRENTbytes = f.read(4)
        if len(BPOPENCURRENTbytes) != 4 :
            print("err:BPOPENCURRENTbytes")
            return "err:BPOPENCURRENTbytes\r\n";
        BPOPENPRESSbytes = f.read(4)
        if len(BPOPENPRESSbytes) != 4 :
            print("err:BPOPENPRESSbytes")
            return "err:BPOPENPRESSbytes\r\n";
        timestamp = 0
        bpopen_volt = 0
        bpopen_current = 0
        bpopen_press = 0
        try :
            # Unpack Integer of 4 bytes
            timestamp = struct.unpack('<I', TIMESTAMPbytes)[0]
            bpopen_volt = struct.unpack('<l', BPOPENVOLTbytes)[0]
            bpopen_current = struct.unpack('<l', BPOPENCURRENTbytes)[0]
            bpopen_press = struct.unpack('<l', BPOPENPRESSbytes)[0]
        except :
            traceback.print_exc()
            print("err:UNPACKVALVEPOWER")
            return "err:UNPACKVALVEPOWER\r\n"
        format = str(timestamp) + ":[BYPASS,0104]battery %5dmV, %7duA (bypass opening), P%+7dmbar\r\n"
        return format % (bpopen_volt,bpopen_current,bpopen_press)
    elif shortId == 7 :
        # Bpclose power LOG
        TIMESTAMPbytes = f.read(4)
        if len(TIMESTAMPbytes) != 4 :
            print("err:TIMESTAMPbytes")
            return "err:TIMESTAMPbytes\r\n";
        BPCLOSEVOLTbytes = f.read(4)
        if len(BPCLOSEVOLTbytes) != 4 :
            print("err:BPCLOSEVOLTbytes")
            return "err:BPCLOSEVOLTbytes\r\n";
        BPCLOSECURRENTbytes = f.read(4)
        if len(BPCLOSECURRENTbytes) != 4 :
            print("err:BPCLOSECURRENTbytes")
            return "err:BPCLOSECURRENTbytes\r\n";
        BPCLOSEPRESSbytes = f.read(4)
        if len(BPCLOSEPRESSbytes) != 4 :
            print("err:BPCLOSEPRESSbytes")
            return "err:BPCLOSEPRESSbytes\r\n";
        timestamp = 0
        bpclose_volt = 0
        bpclose_current = 0
        bpclose_press = 0
        try :
            # Unpack Integer of 4 bytes
            timestamp = struct.unpack('<I', TIMESTAMPbytes)[0]
            bpclose_volt = struct.unpack('<l', BPCLOSEVOLTbytes)[0]
            bpclose_current = struct.unpack('<l', BPCLOSECURRENTbytes)[0]
            bpclose_press = struct.unpack('<l', BPCLOSEPRESSbytes)[0]
        except :
            traceback.print_exc()
            print("err:UNPACKVALVEPOWER")
            return "err:UNPACKVALVEPOWER\r\n"
        format = str(timestamp) + ":[BYPASS,0106]battery %5dmV, %7duA (bypass closing), P%+7dmbar\r\n"
        return format % (bpclose_volt,bpclose_current,bpclose_press)
    elif shortId == 8 :
        # Mermaid detect
        TIMESTAMPbytes = f.read(4)
        if len(TIMESTAMPbytes) != 4 :
            print("err:TIMESTAMPbytes")
            return "err:TIMESTAMPbytes\r\n";
        timestamp = 0
        try :
            # Unpack Integer of 4 bytes
            timestamp = struct.unpack('<I', TIMESTAMPbytes)[0]
        except :
            traceback.print_exc()
            print("err:UNPACKMERMAIDDETECT")
            return "err:UNPACKMERMAIDDETECT\r\n"
        format = str(timestamp) + ":[MRMAID,0027]0dbar, 0degC\r\n"
        return format
    elif shortId == 9 :
        # SBEx1 measures
        TIMESTAMPbytes = f.read(4)
        if len(TIMESTAMPbytes) != 4 :
            print("err:TIMESTAMPbytes")
            return "err:TIMESTAMPbytes\r\n";
        SBEPRESSUREbytes = f.read(4)
        if len(SBEPRESSUREbytes) != 4 :
            print("err:SBEPRESSUREbytes")
            return "err:SBEPRESSUREbytes\r\n";
        SBETEMPbytes = f.read(4)
        if len(SBETEMPbytes) != 4 :
            print("err:SBETEMPbytes")
            return "err:SBETEMPbytes\r\n";
        SBESALbytes = f.read(4)
        if len(SBESALbytes) != 4 :
            print("err:SBESALbytes")
            return "err:SBESALbytes\r\n";
        timestamp = 0
        sbe_press = 0
        sbe_temp = 0
        sbe_sal = 0
        try :
            # Unpack Integer of 4 bytes
            timestamp = struct.unpack('<I', TIMESTAMPbytes)[0]
            sbe_press = struct.unpack('<l', SBEPRESSUREbytes)[0]
            sbe_temp = struct.unpack('<l', SBETEMPbytes)[0]
            sbe_sal = struct.unpack('<l', SBESALbytes)[0]
        except :
            traceback.print_exc()
            print("err:UNPACKSBE")
            return "err:UNPACKSBE\r\n"
        format = str(timestamp) + ":[SBE61 ,0396]P%+7d,T%+7d,S+7%d\r\n"
        return format % (sbe_press,sbe_temp,sbe_sal)
    elif shortId == 10 :
        # Mermaid start
        TIMESTAMPbytes = f.read(4)
        if len(TIMESTAMPbytes) != 4 :
            print("err:TIMESTAMPbytes")
            return "err:TIMESTAMPbytes\r\n";
        timestamp = 0
        try :
            # Unpack Integer of 4 bytes
            timestamp = struct.unpack('<I', TIMESTAMPbytes)[0]
        except :
            traceback.print_exc()
            print("err:UNPACKMERMAIDSTART")
            return "err:UNPACKMERMAIDSTART\r\n"
        format = str(timestamp) + ":[MRMAID,0002]acq started\r\n"
        return format
    elif shortId == 11 :
        # Mermaid stop
        TIMESTAMPbytes = f.read(4)
        if len(TIMESTAMPbytes) != 4 :
            print("err:TIMESTAMPbytes")
            return "err:TIMESTAMPbytes\r\n";
        timestamp = 0
        try :
            # Unpack Integer of 4 bytes
            timestamp = struct.unpack('<I', TIMESTAMPbytes)[0]
        except :
            traceback.print_exc()
            print("err:UNPACKMERMAIDSTOP")
            return "err:UNPACKMERMAIDSTOP\r\n"
        format = str(timestamp) + ":[MRMAID,0003]acq stopped\r\n"
        return format
    elif shortId == 12 :
        # GPS POS
        TIMESTAMPbytes = f.read(4)
        if len(TIMESTAMPbytes) != 4 :
            print("err:TIMESTAMPbytes")
            return "err:TIMESTAMPbytes\r\n";
        LATHEMIbytes = f.read(1)
        if len(LATHEMIbytes) != 1 :
            print("err:LATHEMIbytes")
            return "err:LATHEMIbytes\r\n";
        LATDEGbytes = f.read(1)
        if len(LATDEGbytes) != 1 :
            print("err:LATDEGbytes")
            return "err:LATDEGbytes\r\n";
        LATMINbytes = f.read(1)
        if len(LATMINbytes) != 1 :
            print("err:LATMINbytes")
            return "err:LATMINbytes\r\n";
        LATMMbytes = f.read(2)
        if len(LATMMbytes) != 2 :
            print("err:LATMMbytes")
            return "err:LATMMbytes\r\n";

        LONGHEMIbytes = f.read(1)
        if len(LONGHEMIbytes) != 1 :
            print("err:LONGHEMIbytes")
            return "err:LONGHEMIbytes\r\n";
        LONGDEGbytes = f.read(1)
        if len(LONGDEGbytes) != 1 :
            print("err:LONGDEGbytes")
            return "err:LONGDEGbytes\r\n";
        LONGMINbytes = f.read(1)
        if len(LONGMINbytes) != 1 :
            print("err:LONGMINbytes")
            return "err:LONGMINbytes\r\n";
        LONGMMbytes = f.read(2)
        if len(LONGMMbytes) != 2 :
            print("err:LONGMMbytes")
            return "err:LONGMMbytes\r\n";

        timestamp = 0
        lat_hemi = '?'
        lat_deg = 0
        lat_min = 0
        lat_mm = 0

        long_hemi = '?'
        long_deg = 0
        long_min = 0
        long_mm = 0
        try :
            # Unpack Integer of 4 bytes
            timestamp = struct.unpack('<I', TIMESTAMPbytes)[0]
            lat_hemi = struct.unpack('<B', LATHEMIbytes)[0]
            lat_deg = struct.unpack('<B', LATDEGbytes)[0]
            lat_min = struct.unpack('<B', LATMINbytes)[0]
            lat_mm = struct.unpack('<H', LATMMbytes)[0]

            long_hemi = struct.unpack('<B', LONGHEMIbytes)[0]
            long_deg = struct.unpack('<B', LONGDEGbytes)[0]
            long_min = struct.unpack('<B', LONGMINbytes)[0]
            long_mm = struct.unpack('<H', LONGMMbytes)[0]
        except :
            traceback.print_exc()
            print("err:UNPACKGPSPOS")
            return "err:UNPACKGPSPOS\r\n"
        format = str(timestamp) + ":[SURF  ,0082]%c%02ddeg%02d.%03dmn, %c%03ddeg%02d.%03dmn\r\n"
        return format % (lat_hemi,lat_deg,lat_min,lat_mm,long_hemi,long_deg,long_min,long_mm)
    elif shortId == 13 :
        # GPS DOP
        TIMESTAMPbytes = f.read(4)
        if len(TIMESTAMPbytes) != 4 :
            print("err:TIMESTAMPbytes")
            return "err:TIMESTAMPbytes\r\n";
        HDOPbytes = f.read(1)
        if len(HDOPbytes) != 1 :
            print("err:HDOPbytes")
            return "err:HDOPbytes\r\n";
        mHDOPbytes = f.read(2)
        if len(mHDOPbytes) != 2 :
            print("err:mHDOPbytes")
            return "err:mHDOPbytes\r\n";
        VDOPbytes = f.read(1)
        if len(VDOPbytes) != 1 :
            print("err:VDOPbytes")
            return "err:VDOPbytes\r\n";
        mVDOPbytes = f.read(2)
        if len(mVDOPbytes) != 2 :
            print("err:mVDOPbytes")
            return "err:mVDOPbytes\r\n";
        timestamp = 0
        hdop = 0
        mhdop = 0
        vdop = 0
        mvdop = 0
        try :
            # Unpack Integer of 4 bytes
            timestamp = struct.unpack('<I', TIMESTAMPbytes)[0]
            hdop = struct.unpack('<b', HDOPbytes)[0]
            mhdop = struct.unpack('<h', mHDOPbytes)[0]
            vdop = struct.unpack('<b', VDOPbytes)[0]
            mvdop = struct.unpack('<h', mVDOPbytes)[0]
        except :
            traceback.print_exc()
            print("err:UNPACKGPSDOP")
            return "err:UNPACKGPSDOP\r\n"
        format = str(timestamp) + ":[SURF  ,0084]hdop %d.%03d, vdop %d.%03d\r\n"
        return format % (hdop,mhdop,vdop,mvdop)

    return ""
# Decrypt one file with LOG, WARN,and ERR cards give in arguments
def decrypt_one(path,LOG_card,WARN_card,ERR_card):
    #parse data
    string =""
    with open(path, "rb") as f:
        byte = f.read(1)
        while byte != b'':
            byte = f.read(1)
            if byte != b'#':
                if byte != b'@':
                    continue
                else :
                    string += decrypt_short(f);
            else :
                byte = f.read(1)
                if byte != b'*':
                    continue
                else :
                    string += decrypt_explicit(f,LOG_card,WARN_card,ERR_card);
    return string
# Decrypt all BIN files in a path
def decrypt_all(path):
    # Generate List of BINS file
    files_to_decrypt = glob.glob(path + "*.BIN")
    files_decrypted = list()
    for binary_file in files_to_decrypt :
        # Get version line
        with open(binary_file, "r", errors='replace') as f:
            version = f.readline()
        # Get version
        catch = re.findall("<BDD ([0-9]{3})\.[0-9]{3}\.[0-9]{3}_?V?([0-9]*\.[0-9]+)-?.*>", version)
        if len(catch) > 0:
            # Get database binary_file path
            file_version=catch[-1][1].split(".")
            file_major = '2'
            file_minor = '17'
            if file_version[0] :
                file_major = file_version[0]
            if file_version[1] :
                file_minor = file_version[1]
            file_version = file_major+'.'+file_minor
            model = 0
            if catch[-1][0] == "589" :
                model = 1
            elif catch[-1][0] == "458" :
                model = 2
            database_file = database_get_version(file_version,model)
            if database_file != "" :
                database_file_path = os.path.join(database_path,database_file)
                if os.path.exists(database_file_path):
                    # Read and Parse Database binary_file
                    database =""
                    with open(database_file_path,"r") as f:
                        database = f.read()
                    bin = b''
                    with open(binary_file,"rb") as bin_file:
                        bin = bin_file.read()

                    log_file = binary_file.replace(".BIN",".LOG")
                    binary_file_name = os.path.basename(binary_file)
                    log_file_name = binary_file_name.replace(".BIN",".LOG")

                    print("convert " + binary_file_name + " to " + log_file_name)
                    decrypt_list = json.loads(database)
                    for decrypt_card in decrypt_list:
                        if decrypt_card["TYPE"] == "LOG":
                            log_card = decrypt_card["DECRYPTCARD"]
                        elif decrypt_card["TYPE"] == "WARN":
                            warn_card = decrypt_card["DECRYPTCARD"]
                        elif decrypt_card["TYPE"] == "ERR":
                            err_card = decrypt_card["DECRYPTCARD"]
                    try :
                        result = decrypt_one(binary_file,log_card,warn_card,err_card)
                    except:
                        print(("FORMAT ERROR :" +str(binary_file)))
                    else:
                        with open(log_file,"w") as f:
                            f.write(result)
                        files_decrypted.append(log_file)
                else:
                    print(("No database : " + str(database_file_path)))
    return files_decrypted