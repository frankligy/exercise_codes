<?php

# Ends with a semicolon;
# Every variable starts with $
# implicit type coersion, ===, !==
/*
1. (int) 9.9 will be 9
2. False will be 0 when convert to numeric, but normally will be nothing when convert to string. ie. echo False
3. 'sam' will be converted nothing if numeric
4. "" is smart, escape would work, variable would expand
   '' is not smart, eacape wouldn't work, variable wouldn't expand
*/

# array can be a python list or python dictionary
# PHP is 0-based language   
$stuff = array('hi','hey');     // [0] => 'hi', [1] => 'hey'
$stuff1 = array('name' => 'chuck',
                'course' => 'PHP');
echo $stuff1['course'];
print_r($stuff1);
var_dump($stuff1);   # more verbose, will print the type of value

foreach($stuff1 as $k => $v) {
    echo 'key=',$k,'val=',$v;
}

# useful functions, see my Evernote notes
explode(' ',$string);    # equivalent to python, string.split(' ')
        
# super Global Array
// $_GET, $_POST, $_REQUEST;    REQUEST merge GET and POST

# PHP function default calling by value, if you wanna call by reference, use &$array.

# check if a built-in function is present or not in current version
function_exists();

# import module
require
include

?>