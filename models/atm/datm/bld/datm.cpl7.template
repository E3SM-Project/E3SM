#! /bin/csh -f

if ($NINST_ATM > 1) then
   set inst_counter = 1
   set inst_string = ""
   while ($inst_counter <= $NINST_ATM)
      set inst_string = $inst_counter
      if ($inst_counter <= 999) set inst_string = "0$inst_string"
      if ($inst_counter <=  99) set inst_string = "0$inst_string"
      if ($inst_counter <=   9) set inst_string = "0$inst_string"
      set inst_string = _${inst_string}
      if ( ! -f "$CASEROOT/user_nl_datm${inst_string}" ) then
         cp $CODEROOT/atm/datm/bld/user_nl_datm $CASEROOT/user_nl_datm${inst_string}
      endif
      @ inst_counter = $inst_counter + 1
   end 
else
   if ( ! -f "$CASEROOT/user_nl_datm" ) then
      cp $CODEROOT/atm/datm/bld/user_nl_datm $CASEROOT/user_nl_datm
   endif
endif


