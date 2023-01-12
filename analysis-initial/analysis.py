"""
220 cases
TO REMOVE: 4QQB-11-CAC (#166), 4BS2-1-UGU (#131)

ATTRACT: 57 successes
  in the top 2M (60 % of all natives)
paste ATTRACT-5A-top2M.txt  natives-5A.txt | awk 'NR > 1 && $1/$2 >= 0.6'  | wc -l


"""

import numpy as np

to_remove = [166,131,218,219,220]
#ranked = np.loadtxt("histo-5A-top500k.txt", skiprows=1).astype(int)
nnat = np.loadtxt("natives-5A.txt", skiprows=1).astype(int)
#ranked = ranked[:48]; nnat = nnat[:48] ###
frac = ranked.astype(float) / nnat[:, None]

success = (frac >= 0.6)
to_select = 4

def select_columns(suc, score, best_score, selection, mincol, to_select):
    colsum = suc[:, mincol:].sum(axis=0)
    cols = np.argsort(-colsum)
    optimum = colsum[cols[:to_select]].sum()
    if score + optimum <= best_score:
        return None, None
    old_best_score = best_score
    all_best_selection = None
    for col0 in cols:
        col = col0 + mincol
        next_mincol = col+1
        new_success = suc[:, col]
        next_score = score + new_success.sum()
        next_suc = suc.copy()
        next_suc[new_success, :] = 0
        next_suc[:, col] = 0
        next_selection = selection + [col]
        if to_select == 1:
            curr_best_score = next_score
            curr_selection = next_selection
        else:
            curr_selection, curr_best_score = select_columns(next_suc, next_score, best_score, next_selection, next_mincol, to_select-1)
        if curr_best_score is not None and curr_best_score > best_score:
            all_best_selection = curr_selection
            best_score = curr_best_score
    if all_best_selection is None:
        if to_select > 6:
            print(to_select, best_score)
        return None, None
    else:
        print(to_select, best_score, old_best_score)
        return all_best_selection, best_score

colsum = success.sum(axis=0)
success2 = success[:, np.argsort(-colsum)]
selected, best = select_columns(success2, 0, 0, [], 0, to_select)
selx = np.argsort(-colsum)[selected]