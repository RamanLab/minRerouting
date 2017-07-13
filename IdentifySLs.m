initCobraToolbox;
load('iAF1260');
model=iAF1260; clear iAF1260;

eliList=model.rxns(find(findExcRxns(model)));

fastSL(model,0.01,3,eliList);