def _gen_code(self):
    self._str_frames = str(self._frames)

    # Inject code
    string = "//BEGIN DECLARATION OF THE ROBOT '"+self.name+"'\n\n"
    string += self.object_3d_base.gen_code(self.name)

    string += "const list_links_"+self.name+" = [];\n\n"


    for i in range(len(self.links)):
        string += self.links[i].gen_code(self.name)
        string += "list_links_"+self.name+".push(link_"+str(i)+"_"+self.name+");\n\n"


    string += '''const var_''' + self.name + ''' = new Robot(object3d_'''+self.name+''', list_links_'''+self.name+''',''' + self._str_frames + ''');
    sceneElements.push(var_''' + self.name + ''');
    //USER INPUT GOES HERE
    '''

    return string




