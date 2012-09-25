;;; madx-mode-el -- Major mode for editing MAD-X files

;; Author: Oscar Roberto Blanco Garcia 
;; email : <oscar.roberto.blanco.garcia@cern.ch>
;; Version: 1.0
;; Created: 17.05.2012
;; Keywords: MAD-X major-mode

;; This program is free software; you can redistribute it and/or
;; modify it under the terms of the GNU General Public License as
;; published by the Free Software Foundation; either version 2 of
;; the License, or (at your option) any later version.

;; This program is distributed in the hope that it will be
;; useful, but WITHOUT ANY WARRANTY; without even the implied
;; warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
;; PURPOSE.  See the GNU General Public License for more details.

;; You should have received a copy of the GNU General Public
;; License along with this program; if not, write to the Free
;; Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
;; MA 02111-1307 USA

;;; Commentary:
;; 
;; This mode is modified from an example used in a tutorial about Emacs
;; mode creation. The tutorial can be found here:
;; http://two-wugs.net/emacs/mode-tutorial.html

;; Add this to your .emacs file to load and bind it to files with extension
;; .madx

;;; Code:

(defgroup madx nil
 "Major mode to edit MAD-X files scripts in emacs"
 :group 'languages
)

(defvar madx-mode-hook nil)

;(defvar madx-mode-map
;  (let ((madx-mode-map (make-keymap)))
;    (define-key madx-mode-map "\C-j" 'newline-and-indent)
;    madx-mode-map)
;  "Keymap for MAD-X major mode")

(add-to-list 'auto-mode-alist '("\\.madx\\'" . madx-mode))

(defconst madx-font-lock-keywords-simul
  (list
   ; These define the beginning and end of each MAD-X entity definition
  '("\\<\\(APERTURE\\|CONSTRAINT\\|DYNAP\\|E\\(?:MIT\\|ND\\(?:_?MATCH\\)\\)\\|G\\(?:LOBAL\\|WEIGHT\\)\\|LINE\\|MA\\(?:KETHIN\\|TCH\\)\\|OBSERVE\\|PTC_\\(?:ALIGN\\|CREATE_\\(?:LAYOUT\\|UNIVERSE\\)\\|DUMPMAPS\\|E\\(?:ND\\|PLACEMENT\\)\\|KNOB\\|MO\\(?:MENTS\\|VE_TO_LAYOUT\\)\\|NORMAL\\|OBSERVE\\|PRINTFRAMES\\|READ_ERRORS\\|S\\(?:E\\(?:LECT\\(?:_MOMENT\\)?\\|T\\(?:ACCEL_METHOD\\|DEBUGLEVEL\\|EXACTMIS\\|FRINGE \\|KNOBVALUE\\|RADIATION\\|SWITCH\\|T\\(?:IME\\|OTALPATH\\)\\)\\)\\|TART\\)\\|T\\(?:RACK\\(?:_\\(?:END\\|LINE\\)\\)?\\|WISS\\)\\|VARYKNOB\\)\\|RUN\\|S\\(?:ELECT_PTC_NORMAL\\|IXTRACK\\|TART\\)\\|TWISS\\|VARY\\|WEIGHT\\)\\>" 
 . font-lock-builtin-face)
 )
 "Highlighting expressions for MAD-X mode (simul).")

(defconst madx-font-lock-keywords-programflow
  (list
  '("\\<\\(ELSE\\(?:IF\\)?\\|IF\\|MACRO\\|WHILE\\)\\>"
  . font-lock-keyword-face)
  )
  "Highlighting expressions for MAD-X mode (programflow).")

(defconst madx-font-lock-keywords-controlstm
  (list
  '("\\<\\(ASSIGN\\|C\\(?:ALL\\|O\\(?:GUESS\\|NST\\)\\|REATE\\)\\|D\\(?:ELETE\\|UMPSEQU\\)\\|E\\(?:OPTION\\|SAVE\\|X\\(?:EC\\|IT\\)\\)\\|FILL\\|HELP\\|OPTION\\|PRINT\\|QUIT\\|RE\\(?:A\\(?:D\\(?:\\(?:MY\\)?TABLE\\)\\|L\\)\\|SBEAM\\|TURN\\)\\|S\\(?:AVE\\(?:BETA\\)?\\|E\\(?:LECT\\|T\\(?:VARS\\)?\\)\\|HOW\\|TOP\\|YSTEM\\)\\|T\\(?:ABSTRING\\|ITLE\\)\\|\\(?:US\\|VALU\\|WRIT\\)E\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for MAD-X mode (controlstm).")

(defconst madx-font-lock-keywords-elements
  (list
  '("\\<\\(BEAMBEAM\\|CRABCAVITY\\|D\\(?:IPEDGE\\|RIFT\\)\\|E\\(?:\\(?:COLLIM\\|LSEPAR\\)ATOR\\)\\|H\\(?:\\(?:KICKE\\|MONITO\\)R\\)\\|INSTRUMENT\\|KICKER\\|M\\(?:A\\(?:RKER\\|TRIX\\)\\|ONITOR\\|ULTIPOLE\\)\\|OCTUPOLE\\|QUADRUPOLE\\|R\\(?:BEND\\|COLLIMATOR\\|FCAVITY\\)\\|S\\(?:BEND\\|EXTUPOLE\\|OLENOID\\|ROTATION\\)\\|V\\(?:\\(?:KICKE\\|MONITO\\)R\\)\\|YROTATION\\)\\>"
  . font-lock-type-face)
  )
  "Highlighting expressions for MAD-X mode (elements).")

(defconst madx-font-lock-keywords-beamspec
  (list
  '("\\<\\(\\(?:RES\\)?BEAM\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for MAD-X mode (beamspec).")

(defconst madx-font-lock-keywords-matchingmet
  (list
  '("\\<\\(JACOBIAN\\|LMDIF\\|MIGRAD\\|SIMPLEX\\|USE_MACRO\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for MAD-X mode (matchingmet).")

(defconst madx-font-lock-keywords-orbit_corr
  (list
  '("\\<\\(C\\(?:O\\(?:PTION\\|RRECT \\)\\|SAVE\\)\\|ENDTRACK\\|S\\(?:ETCORR\\|ODD\\|URVEY\\|XF\\(?:READ\\|WRITE\\)\\)\\|T\\(?:\\(?:OUSCHE\\|RAC\\)K\\)\\|USE\\(?:KICK\\|MONITOR\\)\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for MAD-X mode (orbit_corr).")

(defconst madx-font-lock-keywords-plot
  (list
  '("\\<\\(EPRINT\\|PLOT\\|R\\(?:ESPLOT\\|PLOT\\|TRACKSTORE\\|VIEWER\\)\\|SETPLOT\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for MAD-X mode (plot).")

(defconst madx-font-lock-keywords-seqediting
  (list
  '("\\<\\(CYCLE\\|ENDEDIT\\|FLATTEN\\|INSTALL\\|MOVE\\|RE\\(?:FLECT\\|MOVE\\)\\|SEQEDIT\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for MAD-X mode (seqediting).")

(defconst madx-font-lock-keywords-parameters
  (list
  '("\\<\\(A\\(?:DD\\|N\\(?:GLE\\|H[XY]\\)\\|PER\\(?:TYPE\\|_TOL\\)\\|RE[XY]\\|SCALE\\|T\\|UTOPLACEDOWNSTREAM\\)\\|B\\(?:A\\(?:LANCE\\|R[ES]\\)\\|B\\(?:DIR\\|EAT\\|ORBIT\\|SHAPE\\)\\|CURRENT\\|E\\(?:AM_ENVELOPE\\|TA\\(?:0\\|QFX\\)?\\)\\|ISEC\\|UNCHED\\|[VY]\\)\\|C\\(?:A\\(?:LLS\\|VALL\\)\\|ENTRE\\|H\\(?:ARGE\\|ROM\\)\\|L\\(?:ASS\\|EAR\\|\\(?:IS\\|OSED_ORBI\\)T\\)\\|O\\(?:L\\(?:OUR\\|UMN\\)\\|MMAND\\|ND\\|OL\\|R\\|UPLE\\)\\|URRENT\\)\\|D\\(?:AMP\\|E\\(?:BUGLEVEL\\|TUNE\\)\\|ISTORT[12]\\|K\\(?:[NS]R\\|[NS]\\)\\|P\\(?:AR[XY]\\|[HS]I\\)\\|QF\\|THETA\\|UMP\\|[PS]\\)\\|E\\(?:CHO\\|IGN\\|LEMENT\\(?:NAME\\|_BY_ELEMENT\\)?\\|NDSEQUENCE\\|RRORS_OUT\\|VE\\(?:N\\|RYSTEP\\)\\|X\\(?:ACT\\(?:MATCH\\|_MISS?\\)?\\|\\(?:TENSIO\\)?N\\)\\|YN\\|[12TXY]\\)\\|F\\(?:ASTUNE\\|FILE\\|I\\(?:LE\\|NTX?\\)\\|LAG\\|O\\(?:\\(?:N\\|RMA\\)T\\)\\|R\\(?:EQ\\|INGE\\|OM\\)\\|ULL\\|[TXY]\\)\\|G\\(?:AMMA\\|CS\\|N\\(?:FU\\|UF\\)\\)\\|H\\(?:A\\(?:LO\\(?:FILE\\)?\\|ML\\|RMON\\|XIS\\)\\|COEFF[NS]\\|GAP\\|KICK\\|M\\(?:AX\\|IN\\)\\|YSTER\\|[12]\\)\\|I\\(?:CASE\\|N\\(?:DEX\\|FO\\|ITIAL\\(?:_MA\\(?:P_MANUAL\\|TRIX_\\(?:MANUAL\\|TABLE\\)\\)\\)?\\|TER\\(?:POLATE\\|VAL\\)\\)\\)\\|K\\(?:0S\\|1S\\|2S\\|3S\\|BUNCH\\|E\\(?:EPORBIT\\|YWORD\\)\\|ICK\\|NL\\|S[IL]\\|[0-3NS]\\)\\|L\\(?:A\\(?:BEL\\|G\\)\\|E\\(?:NGTH\\|VEL\\)\\|OWER\\|RAD\\|SCALE\\|WIDTH\\|YAPUNOV\\)\\|M\\(?:A\\(?:D8\\|GNET_NAME\\|RKER_PLOT\\|SS\\|TRIX_MANUAL\\|X\\(?:A\\(?:CCELERATION\\|PER\\)\\|_MULT_ORD\\)\\)\\|ETHOD\\|LIST\\|O\\(?:DEL?\\|MENTS\\|N\\(?:ERROR\\|O\\(?:MIAL\\|N\\)\\|SCALE\\)\\)\\|RE[XY]\\|SCAL[XY]\\|ULT\\(?:IP\\(?:\\(?:L\\|OLE_ORDER_RANG\\)E\\)\\|_AUTO_OFF\\)\\)\\|N\\(?:CO\\|EWNAME\\|O\\(?:LINE\\|PRINT\\|RM\\(?:_\\(?:NO\\|OUT\\)\\)?\\|SIXTRACK\\|T\\(?:\\(?:IT\\|SIMP\\)LE\\)\\|VERSION\\|_FATAL_STOP\\)?\\|PART\\|ST\\|TPSA\\)\\|O\\(?:FFSET\\(?:ELEM\\|_DELTAP\\)\\|N\\(?:E\\(?:PASS\\|TABLE\\)\\|LY\\(?:\\(?:ORIENTA\\|POSI\\)TION\\)\\)\\|PT\\|R\\(?:BIT\\|DER\\)\\|VERWRITE\\)\\|P\\(?:A\\(?:R\\(?:AMETRIC\\|ENT\\|TICLE\\)\\|TTERN\\)\\|C\\|ERIOD\\|HI0\\|IPEFILE\\|LA\\(?:CE\\|NEX?\\)\\|O\\(?:LYNOMIAL\\|ST?\\)\\|RINT_A\\(?:LL\\|T_END\\)\\|SI0\\|TC\\(?:_TABLE\\)?\\)\\|QUANT\\(?:ITY\\|UM\\)\\|R\\(?:A\\(?:DI\\(?:AT\\(?:E\\|ION\\(?:_\\(?:ENERGY_LOSS\\|MODEL1\\|QUADr\\)\\)?\\)\\|US\\)\\|N\\(?:DOM\\|GE\\(?:_PLOT\\)?\\)\\)\\|BARC\\|E\\(?:CLOSS\\|F\\(?:ER\\|FRAME\\|POS\\)\\|\\(?:PEA\\|S\\(?:E\\|OU\\|PLI\\)\\)T\\)\\|MATRIX\\|O\\(?:OTNTUPLE\\|W\\)\\|S\\(?:CALE\\|PLIT\\)\\)\\|S\\(?:E\\(?:CTOR\\(?:FILE\\|MAP\\|_\\(?:NMUL_MAX\\|nMUL\\)\\)\\|ED\\|LECTED\\|QUENCE\\)\\|IG[ETXY]\\|L\\(?:ICE\\(?:_MAGNETS\\)?\\|OPE\\)\\|P\\(?:ACE_CHARGE\\|EC\\|LIT\\)\\|SCALE\\|T\\(?:A\\(?:RT_STOP\\|TUS\\)\\|EP\\|RATEGY\\|YLE\\)\\|UMM\\|YM\\(?:BOL\\|P\\(?:L\\|RINT\\)\\)\\)\\|T\\(?:A\\(?:BLE\\(?:ALLSTEPS\\)?\\|RGET\\)\\|E\\(?:LL\\|XT\\)\\|H\\(?:ETA0\\|IN\\(?:_FOC\\)?\\)\\|I\\(?:LT\\|ME\\)\\|O\\(?:L\\(?:ERANCE\\)?\\|TALPATH\\)?\\|R\\(?:\\(?:AC\\(?:KFIL\\)?\\|U\\(?:EPROFIL\\|STRANG\\)\\)E\\)\\|URNS\\|YPE\\)\\|U\\(?:PPER\\|SE\\(?:ORBIT\\|_PTCKNOBS\\)\\)\\|V\\(?:AXIS[1-4]?\\|ERIFY\\|KICK\\|M\\(?:AX\\|IN\\)\\|OLT\\)\\|W\\(?:ARN\\(?:ING\\)?\\|IDTH\\)\\|X\\(?:0\\|BEND\\|DISTR\\|MA\\|SIZE\\)\\|Y\\(?:0\\|DISTR\\|MA\\|SIZE\\)\\|Z\\(?:0\\|DISTR\\)\\|lagf\\|r\\(?:ph[12]\\|v[1-4]\\)\\|[HL]\\)\\>"
  . font-lock-doc-face)
  )
  "Highlighting expressions for MAD-X mode (parameters).")

(defconst madx-font-lock-keywords-errordef
  (list
  '("\\<\\(E\\(?:ALIGN\\|FCOMP\\|RROR\\)\\|SETERR\\)\\>"
  . font-lock-warning-face)
  )
  "Highlighting expressions for MAD-X mode (errordef).")

(defconst madx-font-lock-keywords-constants
  (list
  '("\\<\\(ANTIPROTON\\|C\\(?:ENTRE\\|IRCLE\\|LIGHT\\|OLLIM\\)\\|DEGRAD\\|E\\(?:L\\(?:ECTRON\\|LIPSE\\)\\|MASS\\|NTRY\\|XIT\\)?\\|FALSE\\|LHCSCREEN\\|M\\(?:ARGUERITE\\|UMASS\\)\\|NEGMUON\\|P\\(?:I\\|MASS\\|\\(?:OS\\(?:ITR\\|MU\\)\\|ROT\\)ON\\)\\|QELECT\\|R\\(?:A\\(?:CETRACK\\|DDEG\\)\\|ECT\\(?:\\(?:ANGL\\|ELLIPS\\)E\\)\\)\\|SIMPLE\\|T\\(?:EAPOT\\|RUE\\|WOPI\\)\\)\\>"
  . font-lock-constant-face)
  )
  "Highlighting expressions for MAD-X mode (constants).")

(defconst madx-font-lock-keywords-stringatt
  (list
  '("\\<\\(TITLE\\|system\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for MAD-X mode (stringatt).")

(defconst madx-font-lock-keywords-functions
  (list
  '("\\<\\(A\\(?:BS\\|COS\\|\\(?:SI\\|TA\\)N\\)\\|COSH?\\|EXP\\|GAUSS\\|LOG\\(?:10\\)?\\|RANF\\|S\\(?:INH?\\|QRT\\)\\|T\\(?:ANH?\\|GAUSS\\)\\)\\>"
  . font-lock-function-name-face)
  )
  "Highlighting expressions for MAD-X mode (functions).")

(defconst madx-font-lock-keywords-variables_madx
  (list
  '("\\<\\(A\\(?:LF\\(?:A[1-3][1-3]P?\\|[XY]\\)\\|PER_[1-4]\\)\\|BET\\(?:A[1-3][1-3]P?\\|[XY]\\)\\|D\\(?:D\\(?:P[XY]\\|Q[12]\\|[XY]\\)\\|ELTAP\\|ISP[1-4]\\(?:P[1-3]\\)?\\|MU[XY]\\|P[XY]\\|Q[12]\\|[XY]\\)\\|E\\(?:IGN[1-6][1-6]\\|NERGY\\)\\|GAMA[1-3][1-3]P?\\|K\\(?:\\(?:ICK\\)?[1-6]\\)\\|MU\\(?:[1-3]\\|[XY]\\)\\|N\\(?:1\\(?:[XY]_M\\)?\\|AME\\)\\|ON_\\(?:AP\\|ELEM\\)\\|P\\(?:HI[TXY]?\\|[TXY]N\\|[TXY]\\)\\|Q[12]\\|R\\(?:E\\(?:[1-6][1-6]\\)?\\|M[1-6][1-6]\\|TOL\\|[1-6][1-6]\\)\\|T\\(?:HETA\\|M[1-6][1-6][1-6]\\|N\\|[1-6][1-6][1-6]\\)\\|W[TXY]\\|X\\(?:N\\|TOL\\)\\|Y\\(?:N\\|TOL\\)\\|[RSTXY]\\)\\>"
;  . font-lock-negation-char-face)
  . font-lock-variable-name-face)
  )
  "Highlighting expressions for MAD-X mode (variables_madx).")

(defconst madx-font-lock-special_operators
  (list
   '("\\(->\\|:=\\)"
  . font-lock-warning-face)
  )
  "Highlighting expressions for MAD-X mode (variables_madx).")

(defconst madx-font-lock-special_constants
  (list
   '("\\(#[es]\\)"
  . font-lock-constant-face)
  )
  "Highlighting expressions for MAD-X mode (variables_madx).")


(defconst madx-font-lock-keywords-3
  (append 
     madx-font-lock-special_operators
     madx-font-lock-special_constants
     madx-font-lock-keywords-programflow
     madx-font-lock-keywords-simul
     madx-font-lock-keywords-controlstm
     madx-font-lock-keywords-elements
     madx-font-lock-keywords-beamspec
     madx-font-lock-keywords-matchingmet
     madx-font-lock-keywords-orbit_corr
     madx-font-lock-keywords-plot
     madx-font-lock-keywords-seqediting
     madx-font-lock-keywords-parameters
     madx-font-lock-keywords-errordef
     madx-font-lock-keywords-constants
     madx-font-lock-keywords-stringatt
     madx-font-lock-keywords-functions
     madx-font-lock-keywords-variables_madx
     madx-font-lock-special_operators
     madx-font-lock-special_constants
  )
 "Balls-out highlighting in MAD-X mode.")

(defvar madx-font-lock-keywords madx-font-lock-keywords-3
  "Default highlighting expressions for MAD-X mode.")

(defvar madx-mode-syntax-table
  (let ((madx-mode-syntax-table (make-syntax-table)))
	
    ; This is added so entity names with unde rscores can be more easily parsed
	(modify-syntax-entry ?_ "w" madx-mode-syntax-table)
	(modify-syntax-entry ?. "w" madx-mode-syntax-table)
	
	;  Comment styles are same as C++
	(modify-syntax-entry ?/ ". 124 b" madx-mode-syntax-table)
	(modify-syntax-entry ?* ". 23" madx-mode-syntax-table)
	(modify-syntax-entry ?\n "> b" madx-mode-syntax-table)
	(modify-syntax-entry ?! "< b" madx-mode-syntax-table)
	(modify-syntax-entry ?' "|" madx-mode-syntax-table)
	madx-mode-syntax-table)
  "Syntax table for madx-mode")

;;; ### autoload  
(defun madx-mode ()
  "Major mode for editing MAD-X script files"
  (interactive)
  (kill-all-local-variables)
  (setq mode-name "MAD-X")
  (setq major-mode 'madx-mode)
;  (setq comment-start "!")
;  (use-local-map madx-mode-map)
  (set-syntax-table madx-mode-syntax-table)
  (make-local-variable 'font-lock-defaults)
  (setq font-lock-defaults '(madx-font-lock-keywords nil t))
;; Set up search
  (add-hook 'madx-mode-hook
     (lambda ()  (setq case-fold-search t)))
  (run-hooks 'madx-mode-hook)
)
(provide 'madx-mode)

;;; madx-mode.el ends here

