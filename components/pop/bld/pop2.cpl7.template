#! /bin/csh -f

if ($NINST_OCN > 1) then
   set inst_counter = 1
   set inst_string = ""
   while ($inst_counter <= $NINST_OCN)
      set inst_string = $inst_counter
      if ($inst_counter <= 999) set inst_string = "0$inst_string"
      if ($inst_counter <=  99) set inst_string = "0$inst_string"
      if ($inst_counter <=   9) set inst_string = "0$inst_string"
      set inst_string = _${inst_string}
      if ( ! -f "$CASEROOT/user_nl_pop2${inst_string}" ) then
         cp $CODEROOT/ocn/pop2/bld/user_nl_pop2 $CASEROOT/user_nl_pop2${inst_string}
      endif
      @ inst_counter = $inst_counter + 1
   end 
else
   if ( ! -f "$CASEROOT/user_nl_pop2" ) then
      cp $CODEROOT/ocn/pop2/bld/user_nl_pop2 $CASEROOT/user_nl_pop2
   endif
endif




