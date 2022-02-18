from utils import *

# Set config. Restart animation queue
def _set_ani_frame(self, q=None, htm=None):
    if q is None:
        q = self.q0

    if htm is None:
        htm = self.htm

    n = len(self.links)

    self._frames = []
    self.code = ''
    self.add_ani_frame(0, q, htm)
