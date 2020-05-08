<!DOCTYPE html>
<head><title>Guangyuan Li MD5 Cracker</title></head>
<body>
<h1>MD5 cracker</h1>
<p>This application takes an MD5 hash
of a 4-dighit PIN and 
attempts to hash all two-character combinations
to determine the original 4 digits.</p>
<pre>
Debug Output:
<?php
$goodtext = "Not found";
// If there is no parameter, this code is all skipped
if ( isset($_GET['md5']) ) {
    $time_pre = microtime(true);
    $md5 = $_GET['md5'];

    // This is our alphabet
    $txt = "0123456789";
    $show = 15;

    // Outer loop go go through the alphabet for the
    // first position in our "possible" pre-hash
    // text
    for($i=0; $i<strlen($txt); $i++ ) {
        $ch1 = $txt[$i];   


        for($j=0; $j<strlen($txt); $j++ ) {
            $ch2 = $txt[$j];  

            for($m=0; $m<strlen($txt); $m++ ) {
                $ch3 = $txt[$m];

                for($n=0; $n<strlen($txt); $n++ ) {
                    $ch4 = $txt[$n];

                    $try = $ch1.$ch2.$ch3.$ch4;

                    $check= hash('md5',$try);
                    if ($check == $md5) {
                        $goodtext = $try;
                        break;
                    }

                    if ($show>0) {
                        print "$check $try\n";
                        $show = $show - 1;
                    }
                }
            }


         
        }
    }
    // Compute elapsed time
    $time_post = microtime(true);
    print "Elapsed time: ";
    print $time_post-$time_pre;
    print "\n";
}
?>
</pre>
<!-- Use the very short syntax and call htmlentities() -->
<p>Original Text: <?= htmlentities($goodtext); ?></p>
<form>
<input type="text" name="md5" size="60" />
<input type="submit" value="Crack MD5"/>
</form>
<ul>
<li><a href="crankPassWord.php">Reset</a></li>
</ul>
</body>
</html>