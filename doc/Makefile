# .-----------------------------------------------------------------------.
# | SFcollapse1D                                                          |
# | Gravitational collapse of scalar fields in spherical symmetry         |
# |                                                                       |
# | Copyright (c) 2020, Leonardo Werneck                                  |
# |                                                                       |
# | This program is free software: you can redistribute it and/or modify  |
# | it under the terms of the GNU General Public License as published by  |
# | the Free Software Foundation, either version 3 of the License, or     |
# | (at your option) any later version.                                   |
# |                                                                       |
# | This program is distributed in the hope that it will be useful,       |
# | but WITHOUT ANY WARRANTY; without even the implied warranty of        |
# | MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
# | GNU General Public License for more details.                          |
# |                                                                       |
# | You should have received a copy of the GNU General Public License     |
# | along with this program.  If not, see <https://www.gnu.org/licenses/>.|
# .-----------------------------------------------------------------------.
#
# This is a simple make file that compiles the documentation
all: doc

doc: SFcollapse1D.tex SFcollapse1D.bib
	@echo "Compiling the documentation... \c"
	@cp resources/*.eps .
	@pdflatex SFcollapse1D.tex > /dev/null
	@biber    SFcollapse1D     > /dev/null
	@pdflatex SFcollapse1D.tex > /dev/null
	@pdflatex SFcollapse1D.tex > /dev/null
	@echo "cleaning up... \c"
	@rm -rf *.log *.gz *.aux *.out *.bbl *.bcf *.blg *.xml *.toc *.lof *.lot *-converted-to.pdf resources/*-converted-to.pdf *.eps
	@echo "done!"

clean:
	rm -rf *.log *.gz *.aux *.out *.bbl *.bcf *.blg *.xml *.toc *.lof *.lot *.eps resources/*-converted-to.pdf *.pdf
