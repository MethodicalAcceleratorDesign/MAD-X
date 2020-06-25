;;; madx.el --- Major mode for editing MAD-X files in Emacs

;; Copyright (C) 2013, 2016 Oscar BLANCO
;;               2017, 2018 FSF
;;               2019, 2020 Oscar BLANCO

;; Author     : orblancog <orblancog@gmail.com>
;; Maintainer : orblancog
;; Created    : 18 Nov 2017
;; Keywords   : languages
;; Homepage   : https://github.com/orblancog/mad-x_syntax
;; Version    : 1.9

;; This file is not part of GNU Emacs

;; This program is free software: you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or
;; (at your option) any later version.

;; This program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.

;; You should have received a copy of the GNU General Public License
;; along with this program.  If not, see <https://www.gnu.org/licenses/>.

;;; Commentary:
;; FEATURES in version 1.9
;; * Highlights commands, parameters and special operators in MAD-X 5
;; * If the file extension is '.madx' or '.seq' then the buffer is automatically
;;   highlighted, but any buffer can be highlighted by doing :
;;       `M+X madx-mode` RET
;;   where `M` is the **META** character in Emacs (`M` seems to be **ALT** in Linux)
;;   and RET means press RETURN key
;; * If the limit of 80 characters is desired, they could be highlighted differently
;;   by uncomment the line
;;       ;(require 'whitespace)
;;       ;(setq whitespace-line-column 80) ;; limit line length
;;       ;(setq whitespace-style '(face lines-tail))
;;       ;(add-hook 'madx-mode-hook 'whitespace-mode)
;;   in this (madx.el) file.
;; HOW TO INSTALL IT
;; a) Open this file in emacs and execute
;;      `M+X package-install-file` RET `madx.el` RET
;; b)
;; 1. Copy this file (madx.el) to the highlighting definition folder, e.g.
;;      a)  ~/.emacs.d/lisp/  ---> (Emacs v24.5.X or v25)
;;      b)  ~/.emacs.d/       ---> (Emacs v23.X.X)
;;      c)  ~/                ---> (Emacs v21.X.X)
;;    i.e.
;;      $ cp madx.el ~/.emacs.d/lisp/
;; 2. Edit or create your .emacs file, typically in ~/
;;      adding the following block where the load-path must match point 1.
;;      ;;;; START OF BLOCK TO COPY AND UNCOMMENT
;;      (global-font-lock-mode t);; Enable syntax highlighting
;;      (setq font-lock-maximum-decoration t)
;;      (add-to-list 'load-path "~/.emacs.d");; <--- edit according to 1.
;;      (autoload 'madx-mode "madx" "MADX-mode" t)
;;      (setq auto-mode-alist (append '(("\\.\\(\\(madx\\)\\|\\(seq\\)\\)$" . madx-mode))
;;        auto-mode-alist))
;;      ;;;; END OF BLOCK TO COPY AND UNCOMMENT
;; You should now restart EMACS in order to reload the environment variables.
;; INFO
;; * For mad instructions, visit
;;   mad.web.cern.ch/mad/
;; * Other syntax highlightings could be found inside the
;;   mad sources.  Check the 'syntax' folder in the madx dir !
;;   Write me to the email address above about any bug including an example.

;;; History:
;; v 1.0 First release at CERN. File is also available in the
;;       MAD-X sources "syntax" folder.
;; v 1.1 Adding comments and changing some verbosed names
;; v 1.2 email update oscar.roberto.blanco.garcia@cern.ch deprecated
;;       Adding some variables from MAD-X 5.02.10 manual
;;       Cleaning up faces 8D
;;       when exceeding 80 chars->extra chars in red
;; v 1.3 adding color to numbers and ;
;; v 1.4 Changes to put this file in the emacs elpa repository
;;       Changing email address to orblancog@gmail.com
;; v 1.5 Changing licence
;; v 1.6 Adding automatic syntax highlighting for ".seq"
;; v 1.7 Chaging licence to gpl3 and code-checking following GNU ELPA suggestions
;; v 1.8 Code changes following GNU suggestions
;;       Adding comments about package-install-file
;; v 1.9 Adding few ptc words to the highlighting list

;;; Code:
(defgroup madx nil
  "Major mode to edit MAD-X files in emacs."
  :group 'languages)

(defvar madx-mode-hook nil
  "Hook for madx-model initialize nil.")

;;;; add  80 characters line
;;;; (global-whitespace-mode +1)
;(require 'whitespace)
;(setq whitespace-line-column 80) ;; limit line length
;(setq whitespace-style '(face lines-tail))
;(add-hook 'madx-mode-hook 'whitespace-mode)

(defconst madx-font-lock-keywords-face-all
  ;; madx-font-lock-keywords-programflow
  `((,(regexp-opt '(;; font-lock-keyword-face
		    ;; madx-font-lock-keywords-programflow
		   "IF"
		   "ELSEIF"
		   "ELSE"
		   "WHILE"
		   "MACRO"
		   )
		  'words)
     .  font-lock-keyword-face))
  "Highlighting expressions for MAD-X mode (keywords-all).")

(defconst madx-font-lock-builtin-face-all
  ;; madx-font-lock-keywords-tableaccs
  ;; madx-font-lock-keywords-simul
  ;; madx-font-lock-keywords-controlstm
  ;; madx-font-lock-keywords-filehandstm
  ;; madx-font-lock-keywords-tablehandstm
  ;; madx-font-lock-keywords-beamhandstm
  ;; madx-font-lock-keywords-seqediting
  ;; madx-font-lock-keywords-othrcmds
  ;; madx-font-lock-keywords-matchingmet
  ;; madx-font-lock-keywords-orbit_corr
  ;; madx-font-lock-keywords-plot
  ;; madx-font-lock-keywords-stringatt
  `((,(regexp-opt '(;; font-lock-builtin-face
		   ;;  madx-font-lock-keywords-tableaccs
		   "TABLE"
		   "TABINDEX"
		   "TABSTRING"
		   ;; madx-font-lock-keywords-simul
		   "TWISS"
		   "IBS"
		   "LINE"
		   "MAKETHIN"
		   "APERTURE"
		   "SIXTRACK"
		   "DYNAP"
		   "EMIT"
		   "MATCH"
		   "ENDMATCH"
		   "VARY"
		   "CONSTRAINT"
		   "WEIGHT"
		   "GLOBAL"
		   "GWEIGHT"
		   "PTC_TWISS"
		   "PTC_PRINTPARAMETRIC"
		   "PTC_NORMAL"
		   "SELECT_PTC_NORMAL"
		   "PTC_TRACK"
		   "PTC_TRACK_LINE"
		   "PTC_CREATE_UNIVERSE"
		   "PTC_CREATE_LAYOUT"
		   "PTC_READ_ERRORS"
		   "PTC_MOVE_TO_LAYOUT"
		   "PTC_ALIGN"
		   "PTC_END"
		   "PTC_TRACK_END"
		   "START"
		   "RUN"
		   "PTC_OBSERVE"
		   "OBSERVE"
		   "PTC_START"
		   "PTC_SETSWITCH"
		   "PTC_KNOB"
		   "PTC_SETKNOBVALUE"
		   "MATCH WITHPTCKNOBS"
		   "PTC_PRINTFRAMES"
		   "PTC_SELECT"
		   "PTC_SELECT_MOMENT"
		   "PTC_DUMPMAPS"
		   "PTC_EPLACEMENT"
		   "PTC_VARYKNOB"
		   "END_MATCH"
		   "PTC_MOMENTS"
		   "PTC_SETCAVITIES"
		   "PTC_SETDEBUGLEVEL"
		   "PTC_SETACCEL_METHOD"
		   "PTC_SETEXACTMIS"
		   "PTC_SETRADIATION"
		   "PTC_SETTOTALPATH"
		   "PTC_SETTIME"
		   "PTC_SETFRINGE"
		   ;; madx-font-lock-keywords-controlstm              
		   "EXIT"
		   "QUIT"
		   "STOP"
		   "HELP"
		   "SHOW"
		   "VALUE"
		   "OPTION"
		   "EXEC"
		   "SET"
		   "SYSTEM"
		   "TITLE"
		   "USE"
		   "SELECT"
		   ;;  madx-font-lock-keywords-filehandstm
		   "ASSIGN"
		   "CALL"
		   "RETURN"
		   "PRINT"
		   "PRINTF"
		   "RENAMEFILE"
		   "COPYFILE"
		   "REMOVEFILE"
		   ;;  madx-font-lock-keywords-tablehandstm
		   "CREATE"
		   "DELETE"
		   "READTABLE"
		   "READMYTABLE"
		   "WRITE"
		   "SETVARS"
		   "SETVARS_LIN"
		   "FILL"
		   "SHRINK"
		   ;;  madx-font-lock-keywords-beamhandstm
		   "BEAM"
		   "RESBEAM"
		   ;; madx-font-lock-keywords-seqediting
		   "SEQEDIT"
		   "FLATTEN"
		   "CYCLE"
		   "REFLECT"
		   "INSTALL"
		   "MOVE"
		   "REMOVE"
		   "REPLACE"
		   "EXTRACT"
		   "ENDEDIT"
		   "SAVE"        
		   "DUMPSEQU"
		   ;; madx-font-lock-keywords-othrcmds
		   "SAVEBETA"
		   "COGUESS"
		   "CONST"
		   "EOPTION"
		   "ESAVE"
		   "REAL"
		   ;; madx-font-lock-keywords-matchingmet
		   "LMDIF"
		   "MIGRAD"
		   "SIMPLEX"
		   "JACOBIAN"
		   "USE_MACRO"
		   ;; madx-font-lock-keywords-orbit_corr
		   "CORRECT"
		   "USEMONITOR"
		   "USEKICK"
		   "CSAVE"
		   "SETCORR"
		   "COPTION"
		   "SODD"
		   "SURVEY"
		   "SXFREAD"
		   "SXFWRITE"
		   "TOUSCHEK"
		   "TRACK"
		   "ENDTRACK"
		   ;; madx-font-lock-keywords-plot
		   "PLOT"
		   "RPLOT"
		   "RVIEWER"
		   "RTRACKSTORE"
		   "RESPLOT"
		   "SETPLOT"
		   "EPRINT"
		   ;; madx-font-lock-keywords-stringatt
		   "TITLE"
		   "system"
		   )
		  'words)
     . font-lock-builtin-face))
  "Highlighting expressions for MAD-X mode (builtin-all).")

(defconst madx-font-lock-type-face-all
  ;; madx-font-lock-keywords-elements
  `((,(regexp-opt '(; font-lock-type-face
		   ;; madx-font-lock-keywords-elements
		   "DRIFT"
		   "QUADRUPOLE"
		   "SEXTUPOLE"
		   "OCTUPOLE"
		   "SOLENOID"
		   "CRABCAVITY"
		   "RFCAVITY"
		   "DIPEDGE"
		   "MULTIPOLE"
		   "COLLIMATOR"
		   "ECOLLIMATOR"
		   "RCOLLIMATOR"
		   "YROTATION"
		   "SROTATION"
		   "TRANSLATION"
		   "CHANGEREF"
		   "MARKER"
		   "RBEND"
		   "SBEND"
		   "DIPEDGE"
		   "HKICKER"
		   "VKICKER"
		   "KICKER"
		   "TKICKER"
		   "ELSEPARATOR"
		   "HMONITOR"
		   "VMONITOR"
		   "MONITOR"
		   "INSTRUMENT"
		   "PLACEHOLDER"
		   "BEAMBEAM"
		   "MATRIX"
		   "NLLENS"
		   "RFMULTIPOLE"
		   "ELSEPARATOR"
		   )
		  'words)
     . font-lock-type-face))
  "Highlighting expressions for MAD-X mode (type-all).")

(defconst madx-font-lock-warning-face-all
  ;; madx-font-lock-keywords-errordef
  `((,(regexp-opt '(; font-lock-warning-face
		   ;; madx-font-lock-keywords-errordef
		   "ERROR"
		   "EALIGN"
		   "EFCOMP"
		   "SETERR"
		   )
		  'words)
     . font-lock-warning-face))
  "Highlighting expressions for MAD-X mode (warning-all).")

(defconst madx-font-lock-special_operators
  ;; madx-font-lock-special_operators
  `((,(regexp-opt '(; font-lock-warning-face
		   ;; madx-font-lock-special_operators
		   ":="
		   "->"
		   )
		  t)
     . font-lock-warning-face))
  "Highlighting expressions for MAD-X mode (special-operators).")

(defconst madx-font-lock-constant-face-all
  ;; madx-font-lock-keywords-constants
  `((,(regexp-opt '(; font-lock-constant-face
		   ;; madx-font-lock-keywords-constants
		   "POSITRON"
		   "ELECTRON"
		   "PROTON"
		   "ANTIPROTON"
		   "POSMUON"
		   "NEGMUON"
		   "ION"
		   "PI"
		   "TWOPI"
		   "DEGRAD"
		   "RADDEG"
		   "E"
		   "EMASS"
		   "PMASS"
		   "NMASS"
		   "MUMASS"
		   "CLIGHT"
		   "QELECT"
		   "HBAR"
		   "ERAD"
		   "PRAD"
		   "TRUE"
		   "FALSE"
		   "SIMPLE"
		   "COLLIM"
		   "TEAPOT"
		   "HYBRID"
		   "ENTRY"
		   "CENTRE"
		   "EXIT"
		   "CIRCLE"
		   "RECTANGLE"
		   "ELLIPSE"
		   "LHCSCREEN"
		   "MARGUERITE"
		   "RECTELLIPSE"
		   "RACETRACK"
		   "OCTAGON"
		   "TERMINAL"
		   )
		  'words)
     . font-lock-constant-face))
  "Highlighting expressions for MAD-X mode (constant-all).")

(defconst madx-font-lock-special_constants
  ;; madx-font-lock-special_constants
  `((,(regexp-opt '(; font-lock-constant-face
		   ;; madx-font-lock-special_constants
		   "#e"
		   "#s"
		   )
		  t)
     . font-lock-constant-face))
  "Highlighting expressions for MAD-X mode (special_constants).")

(defconst madx-font-lock-doc-face-all
  ;; madx-font-lock-keywords-parameters
  `((,(regexp-opt '(;; font-lock-doc-face
		   ;; madx-font-lock-keywords-parameters
		   "NOEXPR" ;2016.08
		   "KNLL" ;2016.08
		   "CNLL" ;2016.08
		   "ROOTMACRO"
		   "MOMENT_S"
		   "MOMENT"
		   "SURVEYALL"
		   "DELTA_DEPENDENCY"
		   "DELTAP_DEPENDENCY"
		   "IGNORE_MAP_ORBIT"
		   "RING_PARAMETERS"
		   "CENTER_MAGNETS"
		   "SUMMARY_FILE"
		   "SUMMARY_TABLE"
		   "MAPTABLE"
		   "DELS"
		   "TXI"
		   "TYI"
		   "TLI"
		   "UPDATE"
		   "NCORR"
		   "SNGVAL"
		   "SNGCUT"
		   "CORRLIM"
		   "TWORING"
		   "UNITS"
		   "CORZERO"
		   "BEAM1TAB"
		   "BEAM2TAB"
		   "EXTERN"
		   "NAME_COL"
		   "X_COL"
		   "Y_COL"
		   "TWISSUM"
		   "LSQ"
		   "MICADO"
		   "SVD"
		   "THICK"
		   "VECTOR"
		   "expr"
		   "VLENGTH"
		   "SLOW"
		   "PATCH_ANG"
		   "PATCH_TRANS"
		   "ADD_ANGLE"
		   "ADD_PASS"
		   "NEXT_SEQU"
		   "ZERO_SUPPR"
		   "N_BESSEL"
		   "RIPKEN"
		   "NO_CAVITY_TOTALPATH"
		   "PNL"
		   "PSL"
		   "MAKEDIPEDGE"
		   "TRUNCATE"
		   "APPEND"
		   "ROW1"
		   "ROW1"
		   "PARAM"
		   "SINKICK"
		   "SINPEAK"
		   "SINTUNE"
		   "SINPHASE"
		   "L"
		   "K1"
		   "K1S"
		   "TILT"
		   "K2"
		   "K2S"
		   "K3"
		   "K3S"
		   "KS"
		   "K3S"
		   "KSI"
		   "VOLT"
		   "LAG"
		   "FREQ"
		   "HARMON"
		   "rv1"
		   "rv2"
		   "rv3"
		   "rv4"
		   "rph1"
		   "rph2"
		   "lagf"
		   "PARTICLE"
		   "MASS"
		   "CHARGE"
		   "ENERGY"
		   "PC"
		   "GAMMA"
		   "BETA"
		   "BRHO"
		   "EX"
		   "EXN"
		   "EY"
		   "EYN"
		   "ET"
		   "SIGT"
		   "SIGE"
		   "KBUNCH"
		   "NPART"
		   "PLANE"
		   "BCURRENT"
		   "CURRENT"
		   "BUNCHED"
		   "RADIATE"
		   "NORM_NO"
		   "BV"
		   "SEQUENCE"
		   "ENDSEQUENCE"
		   "REFER"
		   "MAD8"
		   "GNFU"
		   "KICK"
		   "HKICK"
		   "VKICK"
		   "ANGLE"
		   "K0"
		   "K0S"
		   "RESPLIT"
		   "DAMP"
		   "QUANTUM"
		   "RECLOSS"
		   "ELEMENT_BY_ELEMENT"
		   "NORM"
		   "NORM_OUT"
		   "RADIATION_MODEL1"
		   "RADIATION_ENERGY_LOSS"
		   "RADIATION_QUADr"
		   "BEAM_ENVELOPE"
		   "SPACE_CHARGE"
		   "FX"
		   "FY"
		   "FT"
		   "FFILE"
		   "E1"
		   "E2"
		   "FINT"
		   "FINTX"
		   "HGAP"
		   "H1"
		   "H2"
		   "FLAG"
		   "RANGE"
		   "PLACE"
		   "PARENT"
		   "KEYWORD"
		   "FROM"
		   "AT"
		   "REFPOS"
		   "LENGTH"
		   "EXACT_MIS"
		   "CLASS"
		   "PATTERN"
		   "FILE"
		   "FORMAT"
		   "RBARC"
		   "TWISS_PRINT"
		   "THREADER"
		   "THIN_FOC"
		   "NO_FATAL_STOP"
		   "TEXT"
		   "BARE"
		   "SLICE"
		   "THICK"
		   "COMMAND"
		   "TO"
		   "NEWNAME"
		   "BY"
		   "SELECTED"
		   "H"
		   "LRAD"
		   "KNL"
		   "KSL"
		   "TYPE"
		   "SIGX"
		   "SIGY"
		   "XMA"
		   "YMA"
		   "BBSHAPE"
		   "WIDTH"
		   "BBDIR"
		   "ECHO"
		   "ECHOMACRO"
		   "TRACE"
		   "VERIFY"
		   "PERIOD"
		   "SECTORMAP"
		   "SECTORFILE"
		   "KEEPORBIT"
		   "USEORBIT"
		   "COUPLE"
		   "FULL"
		   "COLUMN"
		   "MOMENTS"
		   "PARAMETRIC"
		   "CLEAR"
		   "POS"
		   "POLYNOMIAL"
		   "MONOMIAL"
		   "PARAMETRIC"
		   "QUANTITY"
		   "ROW"
		   "SEED"
		   "ADD"
		   "INFO"
		   "DEBUG"
		   "VERBOSE"
		   "TELL"
		   "RESET"
		   "WARN"
		   "LABEL"
		   "APERTYPE"
		   "HAXIS"
		   "HMIN"
		   "HMAX"
		   "VAXIS"
		   "VAXIS1"
		   "VAXIS2"
		   "VAXIS3"
		   "VAXIS4"
		   "VMIN"
		   "VMAX"
		   "BARS"
		   "SYMBOL"
		   "NOVERSION"
		   "NO_FATAL_ERROR"
		   "NO_FATAL_STOP"
		   "INTERPOLATE"
		   "NOLINE"
		   "NOTITLE"
		   "MARKER_PLOT"
		   "RANGE_PLOT"
		   "MULTIPLE"
		   "PTC"
		   "PTC_TABLE"
		   "TRACKFILE"
		   "CAVALL"
		   "MULT_AUTO_OFF"
		   "MAX_MULT_ORD"
		   "SPLIT"
		   "RADIUS"
		   "WARNING"
		   "STYLE"
		   "COLOUR"
		   "TURNS"
		   "EVERYSTEP"
		   "ONETABLE"
		   "TABLEALLSTEPS"
		   "GCS"
		   "ROOTNTUPLE"
		   "EXTENSION"
		   "FASTUNE"
		   "MAXAPER"
		   "LYAPUNOV"
		   "ORBIT"
		   "TOL"
		   "DS"
		   "DPHI"
		   "DTHETA"
		   "DPSI"
		   "MREX"
		   "MREY"
		   "MSCALX"
		   "MSCALY"
		   "AREX"
		   "AREY"
		   "ORDER"
		   "DKN"
		   "DKS"
		   "DKNR"
		   "DKSR"
		   "HYSTER"
		   "HCOEFFN"
		   "HCOEFFS"
		   "BETA0"
		   "RMATRIX"
		   "STEP"
		   "CHROM"
		   "LOWER"
		   "UPPER"
		   "SLOPE"
		   "OPT"
		   "CALLS"
		   "NO"
		   "XDISTR"
		   "YDISTR"
		   "ZDISTR"
		   "TOLERANCE"
		   "STRATEGY"
		   "REPEAT"
		   "BISEC"
		   "COOL"
		   "BALANCE"
		   "RANDOM"
		   "MODEL"
		   "TARGET"
		   "MODE"
		   "MONERROR"
		   "MONON"
		   "MONSCALE"
		   "PLANEX"
		   "COND"
		   "RESOUT"
		   "CLIST"
		   "MLIST"
		   "STATUS"
		   "POST"
		   "FONT"
		   "LWIDTH"
		   "APER_TOL"
		   "APER_OFFSET"
		   "HALOFILE"
		   "PIPEFILE"
		   "DQF"
		   "BETAQFX"
		   "DP"
		   "DPARX"
		   "DPARY"
		   "COR"
		   "BBEAT"
		   "NCO"
		   "HALO"
		   "INTERVAL"
		   "SPEC"
		   "NOTSIMPLE"
		   "TRUEPROFILE"
		   "OFFSETELEM"
		   "XSIZE"
		   "YSIZE"
		   "ASCALE"
		   "LSCALE"
		   "SSCALE"
		   "RSCALE"
		   "DETUNE"
		   "DISTORT1"
		   "DISTORT2"
		   "START_STOP"
		   "MULTIPOLE_ORDER_RANGE"
		   "NOPRINT"
		   "PRINT_ALL"
		   "PRINT_AT_END"
		   "NOSIXTRACK"
		   "X0"
		   "Y0"
		   "Z0"
		   "THETA0"
		   "PHI0"
		   "PSI0"
		   "SUMM"
		   "CENTRE"
		   "SECTOR_NMUL_MAX"
		   "SECTOR_nMUL"
		   "NTPSA"
		   "SYMPRINT"
		   "TIME"
		   "METHOD"
		   "NST"
		   "EXACT"
		   "OFFSET_DELTAP"
		   "ERRORS_OUT"
		   "ERRORS_IN"
		   "MAGNET_NAME"
		   "RESPLIT"
		   "THIN"
		   "XBEND"
		   "EVEN"
		   "OVERWRITE"
		   "INDEX"
		   "ONEPASS"
		   "DUMP"
		   "DEBUGLEVEL"
		   "LEVEL"
		   "BBORBIT"
		   "SYMPL"
		   "MAXACCELERATION"
		   "EXACT_MISS"
		   "TOTALPATH"
		   "RADIATION"
		   "FRINGE"
		   "ICASE"
		   "CLOSED_ORBIT"
		   "SLICE_MAGNETS"
		   "INITIAL_MATRIX_TABLE"
		   "MATRIX_MANUAL"
		   "INITIAL_MAP_MANUAL"
		   "INITIAL"
		   "ELEMENT"
		   "TRUSTRANGE"
		   "ANHX"
		   "ANHY"
		   "GNUF"
		   "HAML"
		   "EIGN"
		   "INITIAL_MATRIX_MANUAL"
		   "ELEMENTNAME"
		   "KN"
		   "KS"
		   "EXACTMATCH"
		   "ONLYPOSITION"
		   "ONLYORIENTATION"
		   "AUTOPLACEDOWNSTREAM"
		   "REFFRAME"
		   "USE_PTCKNOBS"
		   )
		  'words)
     . font-lock-doc-face))
  "Highlighting expressions for MAD-X mode (doc-all).")

(defconst madx-font-lock-function-name-face-all
  ;; madx-font-lock-keywords-functions
  `((,(regexp-opt '(;;  font-lock-function-name-face
		   ;; madx-font-lock-keywords-functions
		   "SQRT"
		   "LOG"
		   "LOG10"
		   "EXP"
		   "SIN"
		   "COS"
		   "TAN"
		   "ASIN"
		   "ACOS"
		   "ATAN"
		   "SINH"
		   "COSH"
		   "TANH"
		   "SINC"
		   "ABS"
		   "ERF"
		   "ERFC"
		   "FLOOR"
		   "CEIL"
		   "ROUND"
		   "RANF"
		   "GAUSS"
		   "TGAUSS"
		   "FLAT5"
		   "FLAT56"
		   )
		  'words)
     . font-lock-function-name-face))
  "Highlighting expressions for MAD-X mode (name-all)." )

(defconst madx-font-lock-variable-name-face-all
  ;; madx-font-lock-keywords-variables_madx
  `((,(concat (regexp-opt '(;; font-lock-variable-name-face
			    ;; madx-font-lock-keywords-variables_madx
			    "mvar1"
			    "mvar2"
			    "mvar3"
			    "mvar4"
			    "CIRC"
			    "FREQ0"
			    "DTBYDS"
			    "U0"
			    "QS"
			    "ARAD"
			    "PDAMP"
			    "N1MIN"
			    "Z"
			    "PHI"
			    "PSI"
			    "X"
			    "Y"
			    "BETX"
			    "BETY"
			    "NAME"
			    "S"
			    "k0l"
			    "k1l"
			    "k2l"
			    "k3l"
			    "k4l"
			    "K1"
			    "K2"
			    "K3"
			    "K4"
			    "K5"
			    "K6"
			    "K1L"
			    "K2L"
			    "K3L"
			    "K4L"
			    "K5L"
			    "K6L"       
			    "KICK1"
			    "KICK2"
			    "KICK3"
			    "KICK4"
			    "KICK5"
			    "KICK6"
			    "MU1"
			    "MU2"
			    "MU3"
			    "MUX"
			    "MUY"
			    "PX"
			    "PY"
			    "PT"
			    "DELTAP"
			    "XN"
			    "PXN"
			    "WX"
			    "PHI"
			    "THETA"
			    "PHIX"
			    "YN"
			    "PYN"
			    "WY"
			    "PHIY"
			    "TN"
			    "PTN"
			    "WT"
			    "PHIT"
			    "ALFX"
			    "DX"
			    "DPX"
			    "ALFY"
			    "DY"
			    "DPY"
			    "ENERGY"
			    "DMUX"
			    "DDX"
			    "DDPX"
			    "DMUY"
			    "DDY"
			    "DDPY"
			    "Q1"
			    "Q2"
			    "DQ1"
			    "DQ2"
			    "DDQ1"
			    "DDQ2"
			    "N1"
			    "N1X_M"
			    "N1Y_M"
			    "APER_1"
			    "APER_2"
			    "APER_3"
			    "APER_4"
			    "RTOL"
			    "XTOL"
			    "YTOL"
			    "ON_AP"
			    "ON_ELEM"
			    "LENGTH"
			    "ORBIT5"
			    "BETXMAX"
			    "DXMAX"
			    "DXRMS"
			    "XCOMAX"
			    "XRMS"
			    "BETYMAX"
			    "DYMAX"
			    "DYRMS"
			    "YCOMAX"
			    "YCORMS"
			    "SYNCH_1"
			    "SYNCH_2"
			    "SYNCH_3"
			    "SYNCH_4"
			    "SYNCH_5"
			    "DISTANCE"
			    "LYAPUNOV"
			    "LOGDIST"
			    "LOGTURNS"
			    "RE"
					;		    "RE11";REPLACE BY RE[1-6][1-6]
			    "T"
					;		    "T111";; REPLACE BY T[1-6][1-6][1-6]
					;		    "TM111";; REPLACE BY TM[1-6][1-6][1-6]
					;		    "BETA11"; REPLACE BY BETA[1-3][1-3]
					;		    "BETA11P"; REPLACE BY BETA[1-3][1-3]P
			    "ALFA"
					;		    "ALFA11"; REPLACE BY ALFA[1-3][1-3]
					;		    "ALFA11P"; REPLACE BY ALFA[1-3][1-3]P
			    "GAMMATR"
			    "GAMAX"
			    "GAMAY"
					;		    "GAMA11"; REPLACE BY GAMA[1-3][1-3]
					;		    "GAMA11P"; REPLACE BY GAMA[1-3][1-3]P
					;		    "GAMMA11"; REPLACE BY GAMA[1-3][1-3]
					;		    "DISP1; REPLACE BY DISP[1-4]P?[1-3]?"
					;		    "DISP1P1"; REPLACE BY DISP[1-4]P[1-3]
					;		    "EIGN11"; REPLACE BY EIGN[1-6][1-6]
			    "R"
					;		    "R11";; REPLACE BY R[1-6][1-6] AFTER REGEXP-OPT
					;		    "RM11";; REPLACE BY RM[1-6][1-6] AFTER REGEXP-OPT		    
			    )
			  'words)
	      ;; some variables already optimized
	      "\\|\\<RE[1-6][1-6]\\>"
	      "\\|\\<TM?[1-6][1-6][1-6]\\>"
	      "\\|\\<BETA[1-3][1-3]P?\\>"
	      "\\|\\<ALFA[1-3][1-3]P?\\>"
	      "\\|\\<GAMA[1-3][1-3]P?\\>"
	      "\\|\\<DISP[1-4]P?[1-3]?\\>"
	      "\\|\\<EIGN[1-6][1-6]\\>"
	      "\\|\\<RM?[1-6][1-6]\\>")
     . font-lock-variable-name-face))
  "Highlighting expressions for MAD-X mode (variable-name-all).")

(defconst madx-font-lock-intfp-name-face-all
  ;; madx- fonts for integers and floating point numbers
  (list
   '("\\<\\(\\([0-9]+\\.?[0-9]*\\|\\.[0-9]+\\)\\([eE][+-]?\\([0-9]+\\.?[0-9]*\\|[0-9]*\\.[0-9]+\\)\\)?\\)\\>"
     . font-lock-keyword-face))
  "Highlighting expresssions for MAD-X mode (integers and floats).")

(defconst madx-font-lock-keywords-4
  (append
   madx-font-lock-special_constants
   madx-font-lock-special_operators
   madx-font-lock-keywords-face-all
   madx-font-lock-constant-face-all
   madx-font-lock-function-name-face-all
   madx-font-lock-type-face-all
   madx-font-lock-variable-name-face-all
   madx-font-lock-builtin-face-all
   madx-font-lock-warning-face-all
   madx-font-lock-doc-face-all
   madx-font-lock-intfp-name-face-all)
 "Balls-out highlighting in MAD-X mode.")

(defvar madx-font-lock-keywords madx-font-lock-keywords-4
  "Default highlighting expressions for MAD-X mode.")

(defvar madx-mode-syntax-table
  (let ((madx-mode-syntax-table (make-syntax-table)))
	
  ;; This is added so entity names with underscores and dots can be more easily parsed
  (modify-syntax-entry ?_ "w" madx-mode-syntax-table)
  (modify-syntax-entry ?. "w" madx-mode-syntax-table)
	
  ;;  Comment styles are similar to C++
  (modify-syntax-entry ?/ ". 124 b" madx-mode-syntax-table)
  (modify-syntax-entry ?* ". 23" madx-mode-syntax-table)
  (modify-syntax-entry ?\n "> b" madx-mode-syntax-table)
  (modify-syntax-entry ?! "< b" madx-mode-syntax-table)
  (modify-syntax-entry ?' "|" madx-mode-syntax-table)
     madx-mode-syntax-table)
  "Syntax table for `madx-mode'.")

;;;###autoload
(define-derived-mode madx-mode fundamental-mode "madx"
  "Major mode for editing Methodical Accelerator Design X script files."
  (setq font-lock-defaults '(madx-font-lock-keywords nil t))
  ;; Setting up Imenu
  (setq imenu-generic-expression nil)
  (setq imenu-prev-index-position-function nil)
  (setq imenu-extract-index-name-function nil)
  ;;  (imenu-create-index-function)
  ;; Set up search
  (setq case-fold-search t)
  )
;; Highlighting .madx and .seq buffers
;;;###autoload
(add-to-list 'auto-mode-alist '("\\.\\(madx\\|seq\\)\\'" . madx-mode))
(provide 'madx-mode)
;;; madx.el ends here
