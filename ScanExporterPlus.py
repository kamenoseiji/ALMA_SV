#!/usr/bin/env python

import sys
import os;
import urllib2;
from optparse import OptionParser;

class ScanExporter:
    def __init__(self,uid):
        self.uid=uid.replace("___","://").replace("_","/");
        self.uid_=uid.replace("://","___").replace("/","_");
        self.debug=False;
        if not os.path.isdir(self.uid_):
            self.exportMetadata();

    def exportMetadata(self):
        command="asdmExport -m {0}".format(self.uid);
        os.system(command);

    def asdmAsDir(self,asdm):
        return asdm.replace("://","___").replace("/","_");

    def tableFilename(self,table):
        return "{0}/{1}.xml".format(self.uid_,table);
    
    def readTable(self,table):
        lines=open(self.tableFilename(table)).read().split("\n");
        rows=[];read=False;
        for line in lines:
            if "</row>" in line:
                read=False;
            if read:
                tag=line.split("<")[1].split(">")[0];
                value=line.split(">")[1].split("<")[0].strip(" ");
                rows[-1][tag]=value;
            if "<row>" in line:
                read=True;
                rows.append(dict());
        return rows;

    def getScanTree(self):
        scans=self.readTable("Scan");
        subscans=self.readTable("Subscan");
        tree=[];
        for scan in scans:
            tree.append({"scan":scan,"subscans":[subscan for subscan in subscans if subscan["scanNumber"]==scan["scanNumber"]]});
        return tree;

    def getScanNumbers(self):
        tree=self.getScanTree();
        outtree=dict();
        for scan in tree:
            scannumber=scan["scan"]["scanNumber"];
            if scannumber not in outtree:
                outtree[scannumber]=[];
            for subscan in scan["subscans"]:
                outtree[scannumber].append(subscan["subscanNumber"]);
        return outtree;

    def findBinaryUidsForScan(self,scan):
        filename="{0}/{1}.xml".format(self.uid_,"Main");
        lines=open(filename).read().split("\n");
        uids=[];
        scannumber="";
        for line in lines:
            if "<scanNumber>" in line:
                scannumber=line.split(">")[1].split("<")[0];
            if 'entityId="' in line and scannumber==str(scan):
                uid=line.split('entityId="')[1].split('"')[0];
                uids.append(uid);
        return uids;

    def downloadBinaries(self,uids):
        for uid in uids:
            print "Exporting: "+uid;
            url="http://ngas01.sco.alma.cl:7777/RETRIEVE?file_id={0}".format(uid.replace("uid://",""));
            outfile="{0}/ASDMBinary/{1}".format(self.uid_,self.asdmAsDir(uid));
            if not os.path.isfile(outfile):
                command="wget {0} -O {1}".format(url,outfile);
                print command;
                os.system(command);
                #data=urllib2.urlopen(url).read();
                #with open(outfile,"w") as g:
                #    g.write(data);

    def exportScans(self,scans):
        for scan in scans:
            print "Exporting scan "+str(scan);
            uids=self.findBinaryUidsForScan(scan);
            self.downloadBinaries(uids);
            
    def getScanIntents(self):
        rows=self.readTable("Scan");
        return [rowi["scanIntent"] for rowi in rows];
        
    def getScansForIntent(self,intent):
        intents=self.getScanIntents();
        scans=[i+1 for i in range(len(intents)) if intent in intents[i]];
        return scans;            

    def exportByIntent(self,intents):
        scans=[];
        for intent in intents:
            scans+=self.getScansForIntent(intent);
        self.exportScans(scans);




if __name__=="__main__":
    parser=OptionParser();
    parser.add_option("-u", "--uid", dest="uid",help="UID to be processed",default="");
    parser.add_option("-s", "--scans",dest="scans",help="List of scans",default="");
    parser.add_option("-i", "--intent",dest="intent",help="Intent substring",default="");
    (options,args)=parser.parse_args();
    if options.uid=="":
        print "Please specify a uid to be processed with -u";
        sys.exit(1);
    if options.scans=="" and options.intent=="":
        print "Please specify scans to be exported with -s or intent with -i";
    if not options.scans=="":
        uid=options.uid;
        scans=options.scans.split(",");
        se=ScanExporter(uid);
        se.exportScans(scans);
    if not options.intent=="":
        uid=options.uid;
        intent=options.intent;
        se=ScanExporter(uid);
        scans=se.getScansForIntent(intent);
        se.exportScans(scans);
        
    
    
