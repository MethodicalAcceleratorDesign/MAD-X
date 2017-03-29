" Vim syntax file
" Language:    Mad-X
" Maintainer:  Riccardo de Maria <riccardo.de.maria@cern.ch>
" Last Change: 2013 April 30

" Usage:
"
" add in .vimrc the line
" au BufNewFile,BufRead *.madx,*.seq,*.str,*.mask setf madx
"

syntax clear
syntax case ignore
setlocal ignorecase
setlocal iskeyword=a-z,A-Z,48-57,_,.

syntax match madxComment ?//.*$?
syntax match madxComment /!.*$/
syntax region madxComment start=+/\*+ end=+\*/+

syntax match madxPunt2 /;/

syntax match madxOutput /++++++.*/
syntax match madxOutput /+=+=+=.*/
syntax match madxOutput /  + MAD.\{-} +/
syntax match madxOutput /^enter .* module/
syntax match madxOutput /^ GXPLOT-X11.*/
"syntax match madxOutput2 /[ +=]*.*/


"syntax match madxLineReg /\(^\|;\)[^\/!]\_.\{-};/me=e-1 contains=madxPunt,madxNumber,madxFunc,madxType,madxConstant,madxVariables,madxString,madxCommandReg,madxAttribute,madxComment

syntax region madxLineReg start=/[^!]/ end=/;/ contains=madxCommandReg,madxPunt,madxNumber,madxFunc,madxConstant,madxVariables,madxString
"highlight link madxLineReg Todo

"syntax match madxCommandReg contained /\(^\|:[^=]\)\_.\{-}[,;]/me=e-1 contains=,madxLabel,madxCommand,madxElements,madxPunt,madxNumber,madxFunc,madxOutput

syntax region madxCommandReg contained start=/[A-z]/ end=/;/ contains=madxPunt,madxNumber,madxFunc,madxType,madxConstant,madxVariables,madxString,madxCommandReg,madxAttribute,madxComment,madxCommand, madxElements


"highlight link madxCommandReg Todo

syntax match madxPunt contained /+[^+]/me=e-1
syntax match madxPunt contained /=\|(\|)\|:=\|->\|-/
syntax match madxPunt contained /\\/
" YIL13: This match is breaking // and /* .. */ comments
"syntax match madxPunt contained /\*\|\//
syntax match madxPunt contained /:[^=]/me=e-1
syntax match madxPunt contained /,\|;\|{\|}\|\[\|\]/
syntax match madxNumber contained /\<[0-9\.\-][0-9\.\-+ED]*\>/
syntax region madxString contained start=/"/ skip=/\\"/ end=/"/
syntax match madxLabel contained /[A-z0-9 \.]*:[^=]/me=e-1


"Language elements
syntax keyword madxFlowControl contained IF ELSEIF ELSE WHILE MACRO

syntax match madxFunc contained /ATAN/

syntax keyword madxType contained REAL CONST shared

syntax keyword madxConstant contained proton CIRCLE 1 radius  ELLIPSE RECTANGLE LHCSCREEN MARGUERITE RECTELLIPSE RACETRACK POSITRON ELECTRON PROTON ANTIPROTON POSMUON NEGMUON TWOPI

syntax keyword madxElements BEAMBEAM DRIFT ECOLLIMATOR ELSEPARATOR EUROPEAN FOR HKICKER HMONITOR INSTRUMENT KICKER MARKER MATRIX MONITOR MULTIPOLE NUCLEAR OCTUPOLE ORGANIZATION QUADRUPOLE RBEND RCOLLIMATOR RESEARCH RFCAVITY SBEND SEXTUPOLE SOLENOID SROTATION VKICKER VMONITOR YROTATION

syntax keyword madxCommand contained ASSIGN CALL COGUESS CREATE DUMPSEQU EXEC EXIT FILL HELP OPTION PRINT QUIT READTABLE RETURN SAVE SAVEBETA SELECT SET SHOW STOP SYSTEM TITLE USE VALUE WRITE BEAM RESBEAM PLOT RESPLOT SETPLOT SEQEDIT FLATTEN INSTALL MOVE REMOVE CYCLE REFLECT ENDEDIT SEQUENCE ENDSEQUENCE APERTURE SURVEY MAKETHIN

syntax keyword madxCommand contained MATCH ENDMATCH CONSTRAINT VARY LMDIF MIGRAD SIMPLEX WEIGHT JACOBIAN

syntax keyword madxCommand contained TWISS

" Variables
syntax keyword madxVariables contained betx bety dx dy x y px py s mux muy alfx alfy dpx dpy name parent k0l k1l keywordi n1  on_elem keyword

" Attributes
syntax keyword madxAttribute  contained save -echo echo -info info warn clear full vaxis1 vaxis2 format spec

syntax keyword madxAttribute  contained class column default echo file flag label length level pattern period place pmass range sequence slice table text tolerance
syntax keyword madxAttribute  contained BCURRENT BUNCHED BV CHARGE ENERGY ET ETA EX EXN EY EYN GAMMA KBUNCH MASS NPART PARTICLE PC RADIATE SEQUENCE SIGE SIGT

syntax keyword madxAttribute  contained ascale bars colour default file font haxis hmax hmin If interpolate lscale lwidth multiple noline notitle noversion particle post ptc range rscale sscale style symbol table title trackfile vaxis vmax vmin xsize ysize interval

syntax keyword madxAttribute  contained at by class element flag from pattern range sequence start to cor dp

syntax keyword madxAttribute  contained KICK L

"syntax keyword madxAttribute  contained

"syntax match madxAttribute  contained /, *[A-z0-9]\{-} *=/hs=s+1,me=e-1 contains=madxPunt
"syntax match madxAttribute  contained /, *[A-z0-9]\{-} *:=/hs=s+1,me=e-2 contains=madxPunt
"syntax match madxAttribute  contained /, *-*[A-z0-9]\{-} *[,;]/hs=s+1,me=e-1 contains=madxPunt


"element at from angle L sequence KNL L Lrad flag clear particle energy file column echo info warn refer REFPOS range vaxis label place period KICK KSL save bv warn style title haxis VAXIS1 VAXIS2 step upper lower calls tolerance

"defined by the language"
highlight link madxCommand Statement
highlight link madxPunt Statement
highlight link madxPunt2 Statement
highlight link madxFunc Statement
highlight link madxAttribute Type
highlight link madxElements Macro
highlight link madxVariables Special
highlight link madxType     Type

"defined by the user"
highlight link madxComment Comment
highlight link madxLabel String
highlight link madxString String
highlight link madxConstant Constant
highlight link madxNumber Number

highlight link madxOutput Comment
highlight link madxOutput2 Comment

