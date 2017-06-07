#!/home/linhx/bin/luajit

if #arg == 0 then
    print ("<cmd> e.g. zcmd head -3 file.gz")
    os.exit(1)
end

function main()
    for k, v in ipairs(arg) do 
       if v:match(".gz$") then 
          arg[k] = "<(gunzip -c " .. v .. ")"
       end 
    end 
    local cmd = table.concat(arg," ")
--    print(cmd)
    os.execute( "zsh -c \"" .. cmd .. "\"")
end
main()

