source commands.tcl

module PP\#auto -title "PWSCF GUI: module PP.x" -script {

    readfilter  ::pwscf::ppReadFilter

    #
    # Namelist: INPUTPP
    #

    page extract -name "Specify property to calculate" {
	namelist inputpp -name "INPUTPP" {
	    
	    var prefix \
		-label    "Prefix of punch file saved by program PW.X (prefix):" \
		-widget   [list entrybutton "Prefix ..." [list ::pwscf::selectFileRoot $this prefix]] \
		-fmt      %S -validate string \
		-validate string

	    var outdir {
		-label    "Temporary directory where PW.X files resides (outdir):"
		-widget   entrydirselectquote
		-fmt      %S -validate string
		-validate string
	    }
	    var filplot {
		-label    "Output file that will contain the calculated quantity (filplot):"
		-fmt      %S -validate string
		-validate string
	    }
	    var plot_num {
		-label    "What to calculate (plot_num):"
		-widget   radiobox
		-textvalue {
		    "charge density"
		    "total potential (= V_bare + V_H + V_xc)"
		    "local ionic potential"
		    "local density of states at E_fermi" 
		    "local density of electronic entropy"
		    "STM images"
		    "spin polarization (= rho(up) - rho(down))"
		    "|psi|^2"
		    "|psi|^2 (noncollinear case)"
		    "electron localization function (ELF)"
		    "planar average of all |psi|^2"
		    "integrated local density of states (ILDOS)"
		    "the V_bare + V_H potential"
		    "the electric field potential"
		    "the noncolinear magnetization"
		}
		-value { 0 1 2 3 4 5 6 7 7 8 9 10 11 12 13 }
		-fmt %d
	    }
	    var spin_component {
		-label    "Spin component (spin_component):"
		-widget   optionmenu
		-textvalue {
		    "total charge/potential"
		    "spin up charge/potential"
		    "spin down charge/potential"
		    "charge"
		    "absolute value"
		    "x component of the magnetization"
		    "y component of the magnetization"
		    "z component of the magnetization"
		}
		-value { 0 1 2  0  0 1 2 3 }
	    }	
	    
	    separator -label "--- Options for STM images ---"

	    group stm -name "STM" {
		var sample_bias {
		    -label    "For STM: the bias of the sample [in Ryd] in STM images (sample_bias):"
		    -validate fortranreal
		}
		var stm_wfc_matching {
		    -label     "For STM: wave-function matching (stm_wfc_matching):"
		    -widget    radiobox
		    -textvalue {Yes No}
		    -value     {.true. .false.}
		}
		var z {
		    -label    "For STM: height of matching [in celldm(3) units] (z):"
		    -validate fortranreal
		}
		var dz {
		    -label    "For STM: distance of next STM image calculation (dz):"
		    -validate fortranreal
		}
	    }
	    
	    separator -label "--- Options for |psi|^2 ---"

	    group psi2 -name "Psi2" {
		var kpoint {
		    -label    "For |psi^2|: which k-point (kpoint):"
		    -widget    spinint
		    -validate  posint
		    -fmt       %d
		}	
		var kband {
		    -label    "For |psi^2|: which band (kband):"
		    -widget    spinint
		    -validate  posint
		    -fmt       %d
		}
		var lsign {
		    -label    "For |psi^2| & Gamma: save the sign(psi) (lsign):"
		    -widget    radiobox
		    -textvalue {Yes No}
		    -value     {.true. .false.}
		}
	    }
	    
	    separator -label "--- Options for ILDOS ---"

	    group ildos -name "ILDOS" {
		var emin {
		    -label    "For ILDOS: miminum energy [in eV] (emin):"
		    -validate  fortranreal
		}	
		var emax {
		    -label    "FOR ILDOS: maximum energy [in eV] (emax):"
		    -validate  fortranreal
		}	
	    }
	}
    }

    #
    # Namelist: PLOT
    #

    page chdens -name "Specify Plot " {

	namelist plot -name "PLOT" {
	    
	    var nfile {
		-label    "Number of data files (nfile):"
		-widget   spinint
		-validate posint
		-default  1
	    }
	    
	    dimension filepp {
		-label    "Filenames of data files:"
		-start    1
		-end      1
		-widget   entryfileselectquote
		-validate string
		-fmt      %S
	    }
	    
	    dimension weight {
		-label    "Weighting factors:"
		-start    1
		-end      1
		-widget   entry
		-validate fortranreal
		-default  1.0
	    }
	    
	    separator -label "--- Plot info ---"
	    
	    var fileout -label "Name of output file (fileout):" -validate string
	    
	    var iflag {
		-label     "Dimensionality of plot (iflag):"
		-textvalue {
		    "1D plot, spherical average"
		    "1D plot"
		    "2D plot"
		    "3D plot"
		    "2D polar plot"
		}
		-value  { 0 1 2 3 4 }
		-widget optionmenu
	    }
	    
	    var output_format {
		-label     "Format of the output (output_format):"
		-textvalue {
		    "XCRYSDEN's XSF format"
		    "XCRYSDEN's XSF format (whole unit cell)"
		    "format suitable for gnuplot"
		    "format suitable for contour.x"
		    "format suitable for plotrho"
		    "format suitable for gOpenMol"
		    "Gaussian cube-file format"
		}
		-value     { 3 5 0 1 2 4 6 }
		-widget    optionmenu
	    }
	    
	    separator -label "--- Spanning vectors & origin ---"
	    
	    dimension e1 {
		-label    "1st spanning vector:"
		-validate fortranreal
		-start    1
		-end      3
		-pack     left
	    }
	    
	    dimension e2 {
		-label    "2nd spanning vector"
		-validate fortranreal
		-start    1
		-end      3
		-pack     left
	    }
	    
	    dimension e3 {
		-label    "3rd spanning vector"
		-validate fortranreal
		-start    1
		-end      3
		-pack     left
	    }
	    
	    dimension x0 {
		-label    "Origin of the plot"
		-validate fortranreal
		-start    1
		-end      3
		-pack     left
	    }
	    
	    
	    separator -label "--- Number of points in each direction ---" 
	    
	    group nxnynz -name nxnynz {
		packwidgets left
		var nx -label "nx:" -validate posint -widget spinint
		var ny -label "ny:" -validate posint -widget spinint
		var nz -label "nz:" -validate posint -widget spinint
	    }
	    
	    separator -label "--- Polar plot ---"
	    var radius -label "Radius of the sphere (radius):" -validate real
	    
	}
    }
    # ----------------------------------------------------------------------
    # take care of specialties
    # ----------------------------------------------------------------------
    source pp-event.tcl

    # ------------------------------------------------------------------------
    # source the HELP file
    # ------------------------------------------------------------------------
    source pp-help.tcl
}
