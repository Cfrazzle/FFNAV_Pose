function qm = quat_normalize(qm)

%normalize: Note, comment this out if you don't want it!!
nn = (qm(:,1).^2+qm(:,2).^2+qm(:,3).^2+qm(:,4).^2).^0.5;
qm(:,1) = qm(:,1)./nn;
qm(:,2) = qm(:,2)./nn;
qm(:,3) = qm(:,3)./nn;
qm(:,4) = qm(:,4)./nn;

end