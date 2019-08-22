#This file combines parameters from a config file to create as set of command lines to be used to launch a set of fault injection campaigns 
if __name__ == "__main__":
    import argparse
    import json 
    parser = argparse.ArgumentParser(description="Parse input parameters")

    parser.add_argument("--config", help="Config File of parameters that will be parsed to create command line arguments")

    args=parser.parse_args()

    configurationDict=json.load(open(args.config))
    


    #adjust the file names to incorporate the whole path in the dictionary
    encodedPath=configurationDict["encoded_Path"]
    if encodedPath[-1]!="/":
        encodedPath+="/"
    del configurationDict["encoded_Path"]
    numberFiles=0
    archFileList=[]
    for archDict in configurationDict["encoded_files"]:
        for fileIndex, fileName in enumerate(archDict["file_names"]):
            filePath=encodedPath+archDict["file_names"][fileIndex]
            archFileList.append("--arch "+archDict["arch"]+" "+"--fileID "+archDict["file_names"][fileIndex].split(".")[0]+" "+filePath)
            numberFiles+=1

    numberCommands=numberFiles
    #calculate the total number of lines
    for parameter in configurationDict:
        if parameter=="encoded_files":
            continue
        else:
            numberCommands*=len(configurationDict[parameter])

    #generate a set of commands 
    commandList=[]
    for i in range(0,numberCommands):
        commandList.append("")

    for parameter in configurationDict:
        if parameter=="encoded_files":
            continue
        else:
            for i in range(0,numberCommands):
                commandList[i]+="--"+parameter+" "+configurationDict[parameter][i%len(configurationDict[parameter])]+" "

    for i in range(0,numberCommands):
        commandList[i]+=archFileList[i%len(archFileList)]

    with open("intermediate.commands","w+") as writeOut:
        for item in commandList:
            writeOut.write(item+'\n')
        
