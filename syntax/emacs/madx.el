;;; madx-mode -- Major mode for editing MAD-X files in Emacs
;;; FEATURES (v 1.3)
;; * Highlights commands, parameters and special operators in MAD-X 5.2.XX 
;; * If the file extension is '.madx' then the buffer is automatically highlighted,
;;   but any buffer can be highlighted by doing :
;;       `M+X madx-mode`
;;   where `M` is the **META** character in Emacs (`M` seems to be **ALT** in Linux)
;; * If the line is more than 80 characters long, the extra characters are
;;   highlighted differently.
;;   If you dont want this limit, comment/delete the line
;;       (setq whitespace-line-column 80) ;; limit line length
;;   in this file.

;;; LICENCE
;;    Copyright (C) 2016  Oscar BLANCO
;;    This program is free software: you can redistribute it and/or modify
;;    it under the terms of the GNU General Public License as published by
;;    the Free Software Foundation, either version 3 of the License, or
;;    (at your option) any later version.
;;    This program is distributed in the hope that it will be useful,
;;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;;    GNU General Public License for more details.
;;    You should have received a copy of the GNU General Public License
;;    along with this program.  If not, see <http://www.gnu.org/licenses/>.

;;;INSTALLATION
;; 1. Depending on your emacs version, copy this file (madx.el) to the
;;    highlighting definition folder, e.g.
;;      a)  ~/.emacs.d/lisp/  ---> (Emacs v24.5.X)
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
;;      (setq auto-mode-alist (append '(("\\.madx$" . madx-mode))
;;        auto-mode-alist))
;;      ;;;; END OF BLOCK TO COPY AND UNCOMMENT 
;; 3. You should now restart emacs in order to reload the environment variables.

;;;INFO
;; * Author: Oscar Roberto BLANCO GARCIA
;;   email : oscar.blancogarcia@lnf.infn.it
;;   version: 1.3
;;   Created: XX102017 (ddmmyyyy)
;;   Keywords: MAD-X major-mode
;; * New code available at
;;   https://github.com/orblancog/mad-x_syntax.git
;; * For mad instructions, visit
;;   mad.web.cern.ch/mad/
;; * Other syntax highlightings could be found inside the
;;   mad sources. Check the 'syntax' folder in the madx dir !
;;   Write me to the email address above about any bug including an example.
;; * One good example to modify this mode :
;;   http://ergoemacs.org/emacs/elisp_syntax_coloring.html

;;;HISTORY
;; v 1.0 First release at CERN. File is also available in the 
;;       MAD-X sources "syntax" folder.
;; v 1.1 Adding comments and changing some verbosed names
;; v 1.2 email update oscar.roberto.blanco.garcia@cern.ch deprecated
;; v 1.2 Adding some variables from MAD-X 5.02.10 manual
;;       Cleaning up faces 8D
;;       when exceeding 80 chars->extra chars in red
;; v 1.3 adding color to numbers and ;

;;; ... Finally the Code :

(defgroup madx nil
 "Major mode to edit MAD-X files scripts in emacs"
 :group 'languages
)

(defvar madx-mode-hook nil)

(add-to-list 'auto-mode-alist '("\\.madx\\'" . madx-mode))

;; add  80 characters line
;;(global-whitespace-mode +1)
(require 'whitespace)
(setq whitespace-line-column 80) ;; limit line length
(setq whitespace-style '(face lines-tail))
(add-hook 'madx-mode-hook 'whitespace-mode)

(defconst madx-font-lock-keywords-face-all
  ;; madx-font-lock-keywords-programflow
  (list
   '("\\<\\(ELSE\\(?:IF\\)?\\|IF\\|MACRO\\|WHILE\\)\\>"
   .  font-lock-keyword-face)
  )
  "Highlighting expressions for MAD-X mode (keywords-all).")

(defconst madx-font-lock-builtin-face-all
  ;;  madx-font-lock-keywords-tableaccs
  ;; madx-font-lock-keywords-simul
  ;; madx-font-lock-keywords-controlstm
  ;;  madx-font-lock-keywords-filehandstm
  ;;  madx-font-lock-keywords-tablehandstm
  ;;  madx-font-lock-keywords-beamhandstm
  ;; madx-font-lock-keywords-seqediting
  ;; madx-font-lock-keywords-othrcmds
  ;; madx-font-lock-keywords-matchingmet
  ;; madx-font-lock-keywords-orbit_corr
  ;; madx-font-lock-keywords-plot
  ;; madx-font-lock-keywords-stringatt
  (list 
  '("\\<\\(A\\(?:PERTURE\\|SSIGN\\)\\|BEAM\\|C\\(?:ALL\\|O\\(?:GUESS\\|NST\\(?:RAINT\\)?\\|P\\(?:TION\\|YFILE\\)\\|RRECT\\)\\|\\(?:REAT\\|SAV\\|YCL\\)E\\)\\|D\\(?:ELETE\\|UMPSEQU\\|YNAP\\)\\|E\\(?:MIT\\|ND\\(?:EDIT\\|MATCH\\|TRACK\\|_MATCH\\)\\|OPTION\\|PRINT\\|SAVE\\|X\\(?:EC\\|\\(?:I\\|TRAC\\)T\\)\\)\\|F\\(?:ILL\\|LATTEN\\)\\|G\\(?:LOBAL\\|WEIGHT\\)\\|HELP\\|I\\(?:BS\\|NSTALL\\)\\|JACOBIAN\\|L\\(?:INE\\|MDIF\\)\\|M\\(?:A\\(?:KETHIN\\|TCH\\(?: WITHPTCKNOBS\\)?\\)\\|IGRAD\\|OVE\\)\\|O\\(?:BSERVE\\|PTION\\)\\|P\\(?:LOT\\|RINTF?\\|TC_\\(?:ALIGN\\|CREATE_\\(?:LAYOUT\\|UNIVERSE\\)\\|DUMPMAPS\\|E\\(?:ND\\|PLACEMENT\\)\\|KNOB\\|MO\\(?:MENTS\\|VE_TO_LAYOUT\\)\\|NORMAL\\|OBSERVE\\|PRINT\\(?:FRAMES\\|PARAMETRIC\\)\\|READ_ERRORS\\|S\\(?:E\\(?:LECT\\(?:_MOMENT\\)?\\|T\\(?:ACCEL_METHOD\\|CAVITIES\\|DEBUGLEVEL\\|EXACTMIS\\|FRINGE\\|KNOBVALUE\\|RADIATION\\|SWITCH\\|T\\(?:IME\\|OTALPATH\\)\\)\\)\\|TART\\)\\|T\\(?:RACK\\(?:_\\(?:END\\|LINE\\)\\)?\\|WISS\\)\\|VARYKNOB\\)\\)\\|QUIT\\|R\\(?:E\\(?:A\\(?:D\\(?:\\(?:MY\\)?TABLE\\)\\|L\\)\\|FLECT\\|MOVE\\(?:FILE\\)?\\|NAMEFILE\\|PLACE\\|S\\(?:BEAM\\|PLOT\\)\\|TURN\\)\\|PLOT\\|TRACKSTORE\\|UN\\|VIEWER\\)\\|S\\(?:AVE\\(?:BETA\\)?\\|E\\(?:LECT\\(?:_PTC_NORMAL\\)?\\|QEDIT\\|T\\(?:CORR\\|PLOT\\|VARS\\(?:_LIN\\)?\\)?\\)\\|H\\(?:OW\\|RINK\\)\\|I\\(?:MPLEX\\|XTRACK\\)\\|ODD\\|T\\(?:ART\\|OP\\)\\|URVEY\\|XF\\(?:READ\\|WRITE\\)\\|YSTEM\\)\\|T\\(?:AB\\(?:INDEX\\|LE\\|STRING\\)\\|ITLE\\|OUSCHEK\\|RACK\\|WISS\\)\\|USE\\(?:KICK\\|MONITOR\\|_MACRO\\)?\\|VA\\(?:LUE\\|RY\\)\\|W\\(?:EIGHT\\|RITE\\)\\|system\\)\\>"
   . font-lock-builtin-face)
  )
  "Highlighting expressions for MAD-X mode (builtin-all).")

(defconst madx-font-lock-type-face-all
  ;; madx-font-lock-keywords-elements
  (list
   '("\\<\\(BEAMBEAM\\|C\\(?:HANGEREF\\|OLLIMATOR\\|RABCAVITY\\)\\|D\\(?:IPEDGE\\|RIFT\\)\\|E\\(?:\\(?:COLLIM\\|LSEPAR\\)ATOR\\)\\|H\\(?:\\(?:KICKE\\|MONITO\\)R\\)\\|INSTRUMENT\\|KICKER\\|M\\(?:A\\(?:RKER\\|TRIX\\)\\|ONITOR\\|ULTIPOLE\\)\\|NLLENS\\|OCTUPOLE\\|PLACEHOLDER\\|QUADRUPOLE\\|R\\(?:BEND\\|COLLIMATOR\\|F\\(?:CAVITY\\|MULTIPOLE\\)\\)\\|S\\(?:BEND\\|EXTUPOLE\\|OLENOID\\|ROTATION\\)\\|T\\(?:KICKER\\|RANSLATION\\)\\|V\\(?:\\(?:KICKE\\|MONITO\\)R\\)\\|YROTATION\\)\\>"
     . font-lock-type-face)
   )
  "Highlighting expressions for MAD-X mode (type-all).")

(defconst madx-font-lock-warning-face-all
  ;; madx-font-lock-keywords-errordef
  (list
   '("\\<\\(E\\(?:ALIGN\\|FCOMP\\|RROR\\)\\|SETERR\\)\\>"
     . font-lock-warning-face)
   )
  "Highlighting expressions for MAD-X mode (warning-all).")

(defconst madx-font-lock-special_operators
  ;; madx-font-lock-special_operators
  (list
   '("\\(;\\|->\\|:=\\)"
  . font-lock-warning-face)
  )
  "Highlighting expressions for MAD-X mode (special-operators).")

(defconst madx-font-lock-constant-face-all
  ;; madx-font-lock-keywords-constants
  (list
   '("\\<\\(ANTIPROTON\\|C\\(?:ENTRE\\|IRCLE\\|LIGHT\\|OLLIM\\)\\|DEGRAD\\|E\\(?:L\\(?:ECTRON\\|LIPSE\\)\\|MASS\\|NTRY\\|RAD\\|XIT\\)?\\|FALSE\\|H\\(?:BAR\\|YBRID\\)\\|ION\\|LHCSCREEN\\|M\\(?:ARGUERITE\\|UMASS\\)\\|N\\(?:EGMUON\\|MASS\\)\\|OCTAGON\\|P\\(?:I\\|MASS\\|OS\\(?:\\(?:ITR\\|MU\\)ON\\)\\|R\\(?:AD\\|OTON\\)\\)\\|QELECT\\|R\\(?:A\\(?:CETRACK\\|DDEG\\)\\|ECT\\(?:\\(?:ANGL\\|ELLIPS\\)E\\)\\)\\|SIMPLE\\|T\\(?:EAPOT\\|RUE\\|WOPI\\)\\)\\>"
     . font-lock-constant-face)
   )
  "Highlighting expressions for MAD-X mode (constant-all).")

(defconst madx-font-lock-special_constants
    ;; madx-font-lock-special_constants
  (list
   '("\\(#[es]\\)"
  . font-lock-constant-face)
  )
  "Highlighting expressions for MAD-X mode (special_constants).")

(defconst madx-font-lock-doc-face-all
  ;; madx-font-lock-keywords-parameters
  (list
   '("\\<\\(A\\(?:DD\\(?:_\\(?:ANGLE\\|PASS\\)\\)?\\|N\\(?:GLE\\|H[XY]\\)\\|P\\(?:ER\\(?:TYPE\\|_\\(?:OFFSET\\|TOL\\)\\)\\|PEND\\)\\|RE[XY]\\|SCALE\\|T\\|UTOPLACEDOWNSTREAM\\)\\|B\\(?:A\\(?:LANCE\\|R[ES]\\)\\|B\\(?:DIR\\|EAT\\|ORBIT\\|SHAPE\\)\\|CURRENT\\|E\\(?:AM\\(?:1TAB\\|2TAB\\|_ENVELOPE\\)\\|TA\\(?:0\\|QFX\\)?\\)\\|ISEC\\|RHO\\|UNCHED\\|[VY]\\)\\|C\\(?:A\\(?:LLS\\|VALL\\)\\|ENT\\(?:ER_MAGNETS\\|RE\\)\\|H\\(?:ARGE\\|ROM\\)\\|L\\(?:ASS\\|EAR\\|\\(?:IS\\|OSED_ORBI\\)T\\)\\|NLL\\|O\\(?:L\\(?:OUR\\|UMN\\)\\|MMAND\\|ND\\|OL\\|R\\(?:RLIM\\|ZERO\\)?\\|UPLE\\)\\|URRENT\\)\\|D\\(?:AMP\\|E\\(?:BUG\\(?:LEVEL\\)?\\|L\\(?:S\\|TA_DEPENDENCY\\)\\|TUNE\\)\\|ISTORT[12]\\|K\\(?:[NS]R\\|[NS]\\)\\|P\\(?:AR[XY]\\|[HS]I\\)\\|QF\\|THETA\\|UMP\\|[PS]\\)\\|E\\(?:CHO\\(?:MACRO\\)?\\|IGN\\|LEMENT\\(?:NAME\\|_BY_ELEMENT\\)?\\|N\\(?:DSEQUENCE\\|ERGY\\)\\|RRORS_\\(?:IN\\|OUT\\)\\|VE\\(?:N\\|RYSTEP\\)\\|X\\(?:ACT\\(?:MATCH\\|_MISS?\\)?\\|\\(?:TE\\(?:NSIO\\|R\\)\\)?N\\)\\|YN\\|[12TXY]\\)\\|F\\(?:ASTUNE\\|FILE\\|I\\(?:LE\\|NTX?\\)\\|LAG\\|O\\(?:\\(?:N\\|RMA\\)T\\)\\|R\\(?:EQ\\|INGE\\|OM\\)\\|ULL\\|[TXY]\\)\\|G\\(?:AMMA\\|CS\\|N\\(?:FU\\|UF\\)\\)\\|H\\(?:A\\(?:LO\\(?:FILE\\)?\\|ML\\|RMON\\|XIS\\)\\|COEFF[NS]\\|GAP\\|KICK\\|M\\(?:AX\\|IN\\)\\|YSTER\\|[12]\\)\\|I\\(?:CASE\\|GNORE_MAP_ORBIT\\|N\\(?:DEX\\|FO\\|ITIAL\\(?:_MA\\(?:P_MANUAL\\|TRIX_\\(?:MANUAL\\|TABLE\\)\\)\\)?\\|TER\\(?:POLATE\\|VAL\\)\\)\\)\\|K\\(?:0S\\|1S\\|2S\\|3S\\|BUNCH\\|E\\(?:EPORBIT\\|YWORD\\)\\|ICK\\|NLL?\\|S[IL]\\|[0-3NS]\\)\\|L\\(?:A\\(?:BEL\\|G\\)\\|E\\(?:NGTH\\|VEL\\)\\|OWER\\|RAD\\|S\\(?:CALE\\|Q\\)\\|WIDTH\\|YAPUNOV\\)\\|M\\(?:A\\(?:D8\\|GNET_NAME\\|KEDIPEDGE\\|PTABLE\\|RKER_PLOT\\|SS\\|TRIX_MANUAL\\|X\\(?:A\\(?:CCELERATION\\|PER\\)\\|_MULT_ORD\\)\\)\\|ETHOD\\|ICADO\\|LIST\\|O\\(?:DEL?\\|MENT\\(?:_?S\\)?\\|N\\(?:ERROR\\|O\\(?:MIAL\\|N\\)\\|SCALE\\)\\)\\|RE[XY]\\|SCAL[XY]\\|ULT\\(?:IP\\(?:\\(?:L\\|OLE_ORDER_RANG\\)E\\)\\|_AUTO_OFF\\)\\)\\|N\\(?:AME_COL\\|CO\\(?:RR\\)?\\|E\\(?:WNAME\\|XT_SEQU\\)\\|O\\(?:EXPR\\|LINE\\|PRINT\\|RM\\(?:_\\(?:NO\\|OUT\\)\\)?\\|SIXTRACK\\|T\\(?:\\(?:IT\\|SIMP\\)LE\\)\\|VERSION\\|_\\(?:CAVITY_TOTALPATH\\|FATAL_\\(?:ERROR\\|STOP\\)\\)\\)?\\|PART\\|ST\\|TPSA\\|_BESSEL\\)\\|O\\(?:FFSET\\(?:ELEM\\|_DELTAP\\)\\|N\\(?:E\\(?:PASS\\|TABLE\\)\\|LY\\(?:\\(?:ORIENTA\\|POSI\\)TION\\)\\)\\|PT\\|R\\(?:BIT\\|DER\\)\\|VERWRITE\\)\\|P\\(?:A\\(?:R\\(?:AM\\(?:ETRIC\\)?\\|ENT\\|TICLE\\)\\|T\\(?:CH_\\(?:ANG\\|TRANS\\)\\|TERN\\)\\)\\|C\\|ERIOD\\|HI0\\|IPEFILE\\|LA\\(?:CE\\|NEX?\\)\\|NL\\|O\\(?:LYNOMIAL\\|ST?\\)\\|RINT_A\\(?:LL\\|T_END\\)\\|S\\(?:I0\\|L\\)\\|TC\\(?:_TABLE\\)?\\)\\|QUANT\\(?:ITY\\|UM\\)\\|R\\(?:A\\(?:DI\\(?:AT\\(?:E\\|ION\\(?:_\\(?:ENERGY_LOSS\\|MODEL1\\|QUADr\\)\\)?\\)\\|US\\)\\|N\\(?:DOM\\|GE\\(?:_PLOT\\)?\\)\\)\\|BARC\\|E\\(?:CLOSS\\|F\\(?:ER\\|FRAME\\|POS\\)\\|\\(?:PEA\\|S\\(?:E\\|OU\\|PLI\\)\\)T\\)\\|I\\(?:NG_PARAMETERS\\|PKEN\\)\\|MATRIX\\|O\\(?:OT\\(?:MACRO\\|NTUPLE\\)\\|W1?\\)\\|SCALE\\)\\|S\\(?:E\\(?:CTOR\\(?:FILE\\|MAP\\|_\\(?:NMUL_MAX\\|nMUL\\)\\)\\|ED\\|LECTED\\|QUENCE\\)\\|I\\(?:G[ETXY]\\|N\\(?:KICK\\|P\\(?:EAK\\|HASE\\)\\|TUNE\\)\\)\\|L\\(?:ICE\\(?:_MAGNETS\\)?\\|O\\(?:PE\\|W\\)\\)\\|NG\\(?:CUT\\|VAL\\)\\|P\\(?:ACE_CHARGE\\|EC\\|LIT\\)\\|SCALE\\|T\\(?:A\\(?:RT_STOP\\|TUS\\)\\|EP\\|RATEGY\\|YLE\\)\\|U\\(?:MM\\(?:ARY_\\(?:\\(?:FI\\|TAB\\)LE\\)\\)?\\|RVEYALL\\)\\|VD\\|YM\\(?:BOL\\|P\\(?:L\\|RINT\\)\\)\\)\\|T\\(?:A\\(?:BLEALLSTEPS\\|RGET\\)\\|E\\(?:LL\\|XT\\)\\|H\\(?:ETA0\\|I\\(?:CK\\|N\\(?:_FOC\\)?\\)\\|READER\\)\\|I\\(?:LT\\|ME\\)\\|LI\\|O\\(?:L\\(?:ERANCE\\)?\\|TALPATH\\)?\\|R\\(?:\\(?:AC\\(?:KFIL\\)?\\|U\\(?:EPROFIL\\|NCAT\\|STRANG\\)\\)E\\)\\|URNS\\|W\\(?:ISS\\(?:UM\\|_PRINT\\)\\|ORING\\)\\|XI\\|Y\\(?:I\\|PE\\)\\)\\|U\\(?:NITS\\|P\\(?:DATE\\|PER\\)\\|SE\\(?:ORBIT\\|_PTCKNOBS\\)\\)\\|V\\(?:AXIS[1-4]?\\|E\\(?:CTOR\\|R\\(?:BOSE\\|IFY\\)\\)\\|KICK\\|LENGTH\\|M\\(?:AX\\|IN\\)\\|OLT\\)\\|W\\(?:ARN\\(?:ING\\)?\\|IDTH\\)\\|X\\(?:0\\|BEND\\|DISTR\\|MA\\|SIZE\\|_COL\\)\\|Y\\(?:0\\|DISTR\\|MA\\|SIZE\\|_COL\\)\\|Z\\(?:0\\|\\(?:DIST\\|ERO_SUPP\\)R\\)\\|expr\\|lagf\\|r\\(?:ph[12]\\|v[1-4]\\)\\|[HL]\\)\\>"
    . font-lock-doc-face)
   )
  "Highlighting expressions for MAD-X mode (doc-all).")

(defconst madx-font-lock-function-name-face-all
  ;; madx-font-lock-keywords-functions
  (list
   '("\\<\\(A\\(?:BS\\|COS\\|\\(?:SI\\|TA\\)N\\)\\|C\\(?:EIL\\|OSH?\\)\\|E\\(?:RFC?\\|XP\\)\\|FL\\(?:AT56?\\|OOR\\)\\|GAUSS\\|LOG\\(?:10\\)?\\|R\\(?:ANF\\|OUND\\)\\|S\\(?:IN[CH]?\\|QRT\\)\\|T\\(?:ANH?\\|GAUSS\\)\\)\\>"
     . font-lock-function-name-face)
   )
  "Highlighting expressions for MAD-X mode (name-all)" )

(defconst madx-font-lock-variable-name-face-all
  ;; madx-font-lock-keywords-variables_madx
  (list
   '("\\<\\(A\\(?:LF\\(?:A[1-3][1-3]P?\\|[AXY]\\)\\|PER_[1-4]\\|RAD\\)\\|BET\\(?:A[1-3][1-3]P?\\|[XY]MAX\\|[XY]\\)\\|CIRC\\|D\\(?:D\\(?:P[XY]\\|Q[12]\\|[XY]\\)\\|ELTAP\\|\IS\\(?:P\\(?:[1-4]P[1-3]\\|[1-4]\\)\\|TANCE\\)\\|MU[XY]\\|P[XY]\\|Q[12]\\|TBYDS\\|X\\(?:MAX\\|RMS\\)\\|Y\\(?:MAX\\|RMS\\)\\|[XY]\\)\\|E\\(?:IGN[1-6][1-6]\\|NERGY\\)\\\|FREQ0\\|GAM\\(?:A\\(?:[1-3][1-3]P?\\|[XY]\\)\\|MA\\(?:[1-3][1-3]\\|TR\\)\\)\\|K\\(?:1L\\|2L\\|3L\\|4L\\|5L\\|6L\\|ICK[1-6]\\|[1-6]\\)\\|L\\(?:ENGTH\\|OG\\(?:DIST\\\|TURNS\\)\\|YAPUNOV\\)\\|MU[123XY]\\|N\\(?:1\\(?:MIN\\|[XY]_M\\)?\\|AME\\)\\|O\\(?:N_\\(?:AP\\|ELEM\\)\\|RBIT5\\)\\|P\\(?:DAMP\\|HI[TXY]?\\|SI\\|[TXY\]N\\|[TXY]\\)\\|Q[12S]\\|R\\(?:[1-6][1-6]\\|E\\(?:[1-6][1-6]\\)?\\|M[1-6][1-6]\\|TOL\\)\\|SYNCH_[1-5]\\|T\\(?:[1-6][1-6][1-6]\\|HETA\\|M[1-6][1-6][1-6]\\|N\\)\\|U0\\|W[TXY]\\|X\\(?:COMAX\\|N\\|RMS\\\|TOL\\)\\|Y\\(?:CO\\(?:MAX\\|RMS\\)\\|N\\|TOL\\)\\|k\\(?:[0-4]l\\)\\|mvar[1-4]\\|[RSTXYZ]\\)\\>"
     . font-lock-variable-name-face)
   )
  "Highlighting expressions for MAD-X mode (variable-name-all).")

(defconst madx-font-lock-intfp-name-face-all
  ;; madx- fonts for integers and floating point numbers
  (list
   '("\\<\\(\\([0-9]+\\.?[0-9]*\\|\\.[0-9]+\\)\\([eE][+-]?\\([0-9]+\\.?[0-9]*\\|[0-9]*\\.[0-9]+\\)\\)?\\)\\>"
     . font-lock-keyword-face)
   )
  "Highlighting expresssions for MAD-X mode (integers and floats)")

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
   madx-font-lock-intfp-name-face-all
  )
 "Balls-out highlighting in MAD-X mode.")

(defvar madx-font-lock-keywords madx-font-lock-keywords-4
  "Default highlighting expressions for MAD-X mode.")

(defvar madx-mode-syntax-table
  (let ((madx-mode-syntax-table (make-syntax-table)))
	
    ; This is added so entity names with underscores and dots can be more easily parsed
	(modify-syntax-entry ?_ "w" madx-mode-syntax-table)
	(modify-syntax-entry ?. "w" madx-mode-syntax-table)
	
	;  Comment styles are similar to C++
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
