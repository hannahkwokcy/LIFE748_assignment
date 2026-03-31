lrp

* fetch 2gqq
* cealign lrp\_top, 2gqq
* bg\_color white    
* set specular, 0
* color gray80, 2gqq
* color cyan, lrp\_top
* ramp\_new lrp\_Electro\_Scale, lrp\_top\_e\_map, \[-5, 0, 5], \[red, white, blue]
* set surface\_color, lrp\_Electro\_Scale, lrp\_top
* set transparency, 0.3
* select motif\_AF, (lrp\_top and resi 40-65)
* select motif\_EXP, (2gqq and chain C and resi 40-65)
* show sticks, motif\_AF
* show sticks, motif\_EXP
* color yellow, motif\_AF
* color yellow, motif\_EXP



stpA

* fetch 2lrx
* cealign stpA\_top, 2lrx
* bg\_color white
* set specular, 0
* color gray80, 2lrx
* color cyan, stpA\_top
* ramp\_new stpA\_Electro\_Scale, stpA\_top\_e\_map, \[-5, 0, 5], \[red, white, blue]
* set surface\_color, stpA\_Electro\_Scale, stpA\_top
* set transparency, 0.3





argR

* fetch 3v4g
* cealign argR\_top, 3v4g
* bg\_color white
* set specular, 0
* color gray80, 3v4g
* color cyan, argR\_top
* ramp\_new argR\_Electro\_Scale, argR\_top\_e\_map, \[-5, 0, 5], \[red, white, blue]
* set surface\_color, argR\_Electro\_Scale, argR\_top
* set transparency, 0.3
* select argR\_helix, (argR\_top and resi 38-48)
* show sticks, argR\_helix
* color yellow, argR\_helix
* select argR\_wing, (argR\_top and resi 55-65)
* color orange, argR\_wing
* set stick\_radius, 0.15





arsC

fetch 1i9d

* fetch 1i9d
* cealign arsC\_top, 1i9d
* bg\_color white
* set specular, 0
* color gray80, 1i9d
* color cyan, arsC\_top
* select arsc\_triad, (arsC\_top and resi 7+80+82)
* show sticks, arsc\_triad
* color yellow, arsc\_triad
* set stick\_radius, 0.18
* select binding\_loop, (arsC\_top and resi 7-15)
* show cartoon, binding\_loop
* color orange, binding\_loop



nikA

* fetch 1uiv
* cealign nikA\_top, 1uiv
* bg\_color white
* set specular, 0
* color gray80, 1uiv
* color cyan, nikA\_top
* select ni\_loop, (nikA\_top and resi 400-425)
* show sticks, ni\_loop
* set stick\_radius, 0.12
* color orange, ni\_loop
* select ni\_ion, (1uiv and name NI) 
* color green, ni\_ion



