# Update collision object pose (for drawing purpose only)
def _add_col_object(self, sim):
    mth = self.fkm()

    for i in range(len(self.links)):
        for j in range(len(self.links[i].col_objects)):
            sim.add(self.links[i].col_objects[j][0])
