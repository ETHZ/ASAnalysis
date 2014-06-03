#! /usr/bin/python

import os, sys
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

	if not os.path.exists(inputfile) :
		print '[ERROR] Inputfile does not exist!'
		sys.exit(1)

	if os.path.exists(outputfile) :
		print '[warning] Outputfile already exists!'

	input = ROOT.TFile.Open('%s' % (inputfile),'read')
	output = ROOT.TFile.Open('%s' % (outputfile),'recreate')

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
