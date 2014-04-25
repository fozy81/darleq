darleq
======

Demo Phytobenthos classification. This is a demo - has not been tested and is not official tool.

An R function based on WFD method outlined here:

WFD UKTAG document: http://www.wfduk.org/sites/default/files/Media/Environmental%20standards/Annex%202%20Rivers%20Macrophytes%20%26%20Phytobenthos%20DARLEQ.pdf
WFD UKTAG document: http://www.wfduk.org/sites/default/files/Media/Environmental%20standards/Annex%2010%20Lakes%20Macrophytes%20and%20Phytobenthos%20DARLEQ.pdf

Based on tool developed by Steve Juggins. based on work by Kelly, M.G. et al.- see above links for further details.

The directory includes a web server wrapper for the function to create a website using the 'shiny' R package based website
The files breakdwon like this:

* darleqFunc.R - function
* DARLEQ2_TAXON_LIST.csv - list of diatom taxa and related scrores. The darleqFunc.R refers to these when making calculations
* testdata.csv - some test data in correct format for function to use.
* ui.R - interface for website
* server.R - sever side for website



