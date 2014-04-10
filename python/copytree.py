#! /usr/bin/python

import ROOT


def copytree(inputfile, outputfile, tree_name, cut_str):

	print ''
	print '============'
	print '| CopyTree |'
	print '============'
	print ''
	print '         Inputfile:  %s' % (inputfile)
	print '         Outputfile: %s' % (outputfile)
	print ''
		
	input = ROOT.TFile.Open('%s' % (inputfile),'read')
	output = ROOT.TFile.Open('%s' % (outputfile),'recreate')

## NOT NEEDED?	input.cd()
## NOT NEEDED?	obj = ROOT.TObject
## NOT NEEDED?	for key in ROOT.gDirectory.GetListOfKeys():
## NOT NEEDED?		input.cd()
## NOT NEEDED?		obj = key.ReadObj()
## NOT NEEDED?		#print obj.GetName()
## NOT NEEDED?		if obj.GetName() == 'tree':
## NOT NEEDED?			continue
## NOT NEEDED?		output.cd()
## NOT NEEDED?		#print key.GetName()
## NOT NEEDED?		obj.Write(key.GetName())

	inputTree = input.Get('%s' % (tree_name))
	nEntries = inputTree.GetEntries()
	output.cd()
	print '[status] copy tree %s with cuts: %s' % (tree_name, cut_str)
	outputTree = inputTree.CopyTree(cut_str)
	kEntries = outputTree.GetEntries()
	print '         before cuts: %d' % (nEntries)
	print '         survived:    %d' % (kEntries)
	outputTree.AutoSave()
#	output.ls()
	print '[status] writing and closing files..'
	output.Write()
	output.Close()
	input.Close()
	print '[status] exit CopyTree\n'
