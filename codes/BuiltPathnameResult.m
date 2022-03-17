function pathname = BuiltPathnameResult(test, foldername, filename)

foldername = convertStringsToChars(foldername);
filename = convertStringsToChars(filename);

if test
    pathname = fullfile(GetCodesMatthieuPath,'results',foldername,'Test',filename);
else
    pathname = fullfile(GetCodesMatthieuPath,'results',foldername,filename);
end
