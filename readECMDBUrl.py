# import urllib

# link = "http://ecmdb.ca/compounds/M2MDB000084#concentrations"
# f = urllib.urlopen(link)
# myfile = f.read()
# print len(myfile)

# print myfile[myfile.find("Concentrations")-120:myfile.find("Concentrations")+120]

import mechanize
br = mechanize.Browser()
url = "http://ecmdb.ca/compounds/M2MDB000084"
response = br.open(url)
print response.read()      # the text of the page

# allow everything to be written to
# br.set_all_readonly(False)    
# ignore robots
# br.set_handle_robots(False)   
# can sometimes hang without this
# br.set_handle_refresh(False)  


# br.addheaders = [('User-agent', 'Firefox')]