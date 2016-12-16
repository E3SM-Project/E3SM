    function applyFilter(filter_text) {

	// applying a filter hides all standard names not matching filter_text
	// if filter_text contains no spaces, it is treated as a regexp
        // otherwise, all substrings must occur somewhere

	var is_match = false;
	var search_type = 'regexp';
	var search_help_text = false;
	var num_matches = 0;
	var is_boolean_and = true;

	search_help_text = (document.getElementById('search_help_text').checked);
	is_boolean_and = (document.getElementById('logical_operator_and').checked);

	if (filter_text.indexOf(' ') == -1) {
	    search_type = 'regexp';
	    var re = new RegExp(filter_text, 'i')
	}
	else {
	    search_type = 'string';
	    var string_parts = filter_text.split(' ');
	}

	allTRs = document.getElementsByTagName('tr');

	for (var i = 0; i < allTRs.length; i++) {
	    curTR = allTRs[i];

	    if (curTR.id != '') {

		if (search_type == 'regexp') {

		    is_match = curTR.id.substring(0, curTR.id.length - 3).match(re);

		    if (search_help_text) {

			var helpText = document.getElementById(curTR.id.substring(0,curTR.id.length - 3) + '_help').innerHTML;
			is_match = is_match || helpText.match(re);
		    }
		}
		else {

		    if (is_boolean_and) {
			var is_name_match = true;
			for (var j = 0; j < string_parts.length && is_name_match; j++) {

			    if (!curTR.id.match(new RegExp(string_parts[j], 'i'))) {
				is_name_match = false;
			    }
			}
		    }
		    else {

			var is_name_match = false;
			for (var j = 0; j < string_parts.length && !is_name_match; j++) {

			    if (curTR.id.substring(0, curTR.id.length - 3).match(new RegExp(string_parts[j], 'i'))) {
				is_name_match = true;
			    }
			}
		    }

		    is_match = is_name_match;

		    if (search_help_text) {
			var helpText = document.getElementById(curTR.id.substring(0,curTR.id.length - 3) + '_help').innerHTML;

			if (is_boolean_and) {
			    var is_help_match = true;

			    for (var j = 0; j < string_parts.length && is_help_match; j++) {

				if (!helpText.match(new RegExp(string_parts[j], 'i'))) {
				    is_help_match = false;
				}
			    }
			}
			else {

			    var is_help_match = false;

			    for (var j = 0; j < string_parts.length && !is_help_match; j++) {

				if (helpText.match(new RegExp(string_parts[j], 'i'))) {
				    is_help_match = true;
				}
			    }
			}

			is_match = is_match || is_help_match;

		    }
		}

		if (!is_match) {
		    curTR.style.display = 'none';
		}
		else {
		    num_matches++;
		    curTR.style.display = '';
		    if (search_help_text) {
			showHelp(curTR.id.substring(0,curTR.id.length - 3));
		    }
		    else {
			hideHelp(curTR.id.substring(0,curTR.id.length - 3));
		    }
		}
	    }
	}

	var filter_matches = document.getElementById('filter_matches');
	var filter_matches_num = document.getElementById('filter_matches_num');
	var filter_matches_query = document.getElementById('filter_matches_query');

	if (filter_text != '') {
	    filter_matches.style.visibility = 'visible';
	    filter_matches_num.innerHTML = num_matches;
	    filter_matches_query.innerHTML = filter_text;
	}
	else {
	    filter_matches.style.visibility = 'hidden';
	}

    } // end function applyFilter()

    function clearFilter() {

	allTRs = document.getElementsByTagName('tr');

	for (var i = 0; i < allTRs.length; i++) {
	    curTR = allTRs[i];
	    if (curTR.id != '') {
		curTR.style.display = '';
		hideHelp(curTR.id.substring(0,curTR.id.length - 3));

	    }
	}

	var filter_matches = document.getElementById('filter_matches');
	filter_matches.style.visibility = 'hidden';

	document.getElementById('filter_text').value = '';
    }

    function expandAllHelp() {

	// expand help text for all visible elements

	allTRs = document.getElementsByTagName('tr');

	for (var i = 0; i < allTRs.length; i++) {
	    curTR = allTRs[i];
	    if (curTR.id != '' && curTR.style.display != 'none') {
		showHelp(curTR.id.substring(0,curTR.id.length - 3));

	    }
	}
    }

    function collapseAllHelp() {

	// collapse help text for all visible elements

	allTRs = document.getElementsByTagName('tr');

	for (var i = 0; i < allTRs.length; i++) {
	    curTR = allTRs[i];
	    if (curTR.id != '' && curTR.style.display != 'none') {
		hideHelp(curTR.id.substring(0,curTR.id.length - 3));

	    }
	}
    }

    function toggleHelp(standard_name) {

	// check for the existence of the help "tr" object for this standard_name

        var helpDiv = document.getElementById(standard_name + '_help');

	if (helpDiv) {

	    if (helpDiv.style.display != 'none') {

		helpDiv.style.display = 'none';

		curArrow = document.getElementById(standard_name + '_arrow');
		curArrow.src = "./images/arrow_right.gif";
	    }
	    else {
		helpDiv.style.display = '';

		curArrow = document.getElementById(standard_name + '_arrow');
		curArrow.src = "./images/arrow_down.gif";
	    }
	}
    }


    function showHelp(standard_name) {

	var helpDiv = document.getElementById(standard_name + '_help');

	if (helpDiv) {

	    helpDiv.style.display = '';
	    curArrow = document.getElementById(standard_name + '_arrow');
	    curArrow.src = "./images/arrow_down.gif";
	}
    }

    function hideHelp(standard_name) {

	var helpDiv = document.getElementById(standard_name + '_help');

	if (helpDiv) {
	    helpDiv.style.display = 'none';
	    curArrow = document.getElementById(standard_name + '_arrow');
	    curArrow.src = "./images/arrow_right.gif";
	}
    }

