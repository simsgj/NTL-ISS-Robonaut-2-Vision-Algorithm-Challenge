<?php
$config = include("bundle_config.php");
$files = $config['files'];
$exts = $config['extensions'];

$code = "";

$declaration_headers = array();

foreach ($exts as $ext) {
    foreach ($files as $file) {
        $fileName = $file.$ext;
        if (file_exists($fileName)) {
            $str = file_get_contents($fileName);
            $strs = explode("\n", $str);
            for($i = 0; $i < count($strs); $i++) {
                if (substr(trim($strs[$i]), 0, 8) == "#include") {
                    if (substr(trim($strs[$i]), 0, 10) == "#include <") {
                        $declaration_headers[trim($strs[$i])] = 1;
                    }
                    
                    array_splice($strs, $i, 1);
                    $i--;
                }
            }
            $str = implode("\n", $strs);
            $code .= $str."\n\n";
        }
    }
}

$code = implode("\n", array_keys($declaration_headers))."\n\n".$code;
  
file_put_contents("bundle.cpp", $code);
echo "Bundle complete !\n";
?>
