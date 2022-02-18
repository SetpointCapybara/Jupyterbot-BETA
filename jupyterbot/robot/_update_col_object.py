# Update collision object pose (for drawing purpose only)
def _update_col_object(self):
    htm = self.fkm(axis="dh")

    for i in range(len(self.links)):
        for j in range(len(self.links[i].col_objects)):
            self.links[i].col_objects[j][0].add_ani_frame(htm[i] @ self.links[i].col_objects[j][1])
